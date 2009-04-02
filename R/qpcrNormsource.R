setClass("qpcrBatch", representation(geneNames = "character", plateIndex = "character", exprs = "matrix", normalized="logical", normGenes="character"))

setMethod("normalize", "qpcrBatch", 
          function(object, method = c("housekeepinggenes", "quantile", "rankinvariant"), ...) {
              method <- match.arg(method)
              switch(method, 
                     housekeepinggenes = normQpcrHouseKeepingGenes(object, ...),
                     quantile = normQpcrQuantile(object, ...),
                     rankinvariant = normQpcrRankInvariant(object, ...))
          }
)

`calcCV` <-
function(qBatch){
	m <- mean(apply(qBatch@exprs, 1, function(x){ sd(x, na.rm=T)/mean(x, na.rm=T) }), na.rm=T)
	return(m)
}

`ctQc` <-
function(x){ 
	qc.x <- NULL 
	for( i in 1:nrow(x) ){
		ct.trip <- sort(x[i,])
		diff <- double(2) ; diff.reg <- double(2)
		diff[1] <- ct.trip[2] - ct.trip[1] ; diff[2] <- ct.trip[3] - ct.trip[2]
		
		for( j in 1:length(diff) ){
			if( diff[j] <= 0.2 ){ diff.reg[j] <- 1 }
			else if( diff[j] > 0.2 & diff[j] <= 1 ){ diff.reg[j] <- 2 }
			else if( diff[j] > 1 ){ diff.reg[j] <- 3 }
			else{ diff.reg[j] <- 4 }
		}
		
		if( diff[1] == diff[2] ){ ct.fin <- mean(ct.trip, na.rm=TRUE) }
		else if( diff[1] < diff[2] ){ ct.fin <- mean(ct.trip[1:2], na.rm=TRUE) }
		else{ ct.fin <- mean(ct.trip[2:3], na.rm=TRUE) }
	
		qc.x <- c(qc.x, ct.fin) 
	}	
	return(qc.x)
}

`matrixByPlate` <-
function(xvec, plateIndex){

	x <- NULL ; numPlates <- length(unique(plateIndex)) ; uniPlate <- unique(plateIndex) 
	for( i in 1:numPlates ){ 
		x <- c(x, length(plateIndex[plateIndex == uniPlate[i]]))
	}
	
	numDiffSizes <- length(table(x)) 
	matRow <- max(x) ; matCol <- numPlates
	res.mat <- NULL ; res.mat.ind <- NULL 
	
	for( i in 1:numPlates ){
		numGenes <- length(plateIndex[plateIndex == uniPlate[i]])
		xpad <- c(xvec[plateIndex == uniPlate[i]], rep(NA, matRow - numGenes))
		res.mat <- cbind(res.mat, xpad) 
		res.mat.ind <- cbind(res.mat.ind, c(rep(TRUE, numGenes), rep(FALSE, matRow - numGenes)))
	}
	return(list(resMat = res.mat, resMatIndex = res.mat.ind))
	# resMatIndex is matrix of T/F flags where T = not na
}


`normQpcrHouseKeepingGenes` <-
function(qBatch, hkeep.genes){
	avg.dat <- qBatch@exprs 
	norm.xdat <- NULL
	plate.ind <- qBatch@plateIndex 
	if( sum(hkeep.genes %in% qBatch@geneNames) == 0 ){
		notfound <- hkeep.genes[!hkeep.genes %in% qBatch@geneNames]
		if( length(notfound) == 1 ){	# if just 1 gene not found
			e <- simpleError(paste(notfound, " not found in geneNames slot", sep=""))
		}
		else{	# if a couple of genes, insightful to say which ones weren't found
			e <- simpleError(paste(paste(notfound, collapse=", "), " not found in geneNames slot", sep=""))
		}
		stop(e)
	}
	else{
		if( length(hkeep.genes) > 1 | sum(qBatch@geneNames %in% hkeep.genes) > 1 ){
			hkeep.expr <- apply(avg.dat[qBatch@geneNames %in% hkeep.genes,], 2, mean, na.rm=TRUE)
		}
		else{	
			hkeep.expr <- avg.dat[qBatch@geneNames %in% hkeep.genes,]
		}
		scale.hkeep <- hkeep.expr/hkeep.expr[1]
		norm.dat <- avg.dat 
		for( i in 1:ncol(avg.dat) ){
			norm.dat[,i] <- avg.dat[,i]*scale.hkeep[i]
		}
		colnames(norm.dat) <- colnames(avg.dat) ; rownames(norm.dat) <- qBatch@geneNames 
		qBatch@exprs <- norm.dat
		qBatch@normalized <- TRUE 
		qBatch@normGenes <- as.character(hkeep.genes) 
		return(qBatch)
	}
}


`normQpcrQuantile` <-
function(qBatch){
	require(limma)
	avg.dat <- qBatch@exprs	; plateIndex <- qBatch@plateIndex 
	numPlates <- length(unique(plateIndex))
	plateInfo <- table(as.numeric(table(plateIndex)))
	numDiffSizes <- length(plateInfo)
	all.samp.norm <- NULL 
	
	for( i in 1:ncol(avg.dat) ){
		xmat <- matrixByPlate(avg.dat[,i], plateIndex)
		xnorm <- normalizeBetweenArrays(xmat$resMat, method="quantile")
		xvec <- NULL 
		for( j in 1:ncol(xnorm) ){ 
			xvec <- c(xvec, xnorm[,j][xmat$resMatIndex[,j]])
		}
		all.samp.norm <- cbind(all.samp.norm, xvec)
	}
	all.samp.norm <- normalizeBetweenArrays(all.samp.norm, method="quantile")
	rownames(all.samp.norm) <- qBatch@geneNames 
	colnames(all.samp.norm) <- colnames(avg.dat) 
	qBatch@exprs <- all.samp.norm
	qBatch@normalized <- TRUE
	qBatch@normGenes <- qBatch@geneNames
	return(qBatch)
}

`normQpcrRankInvariant` <-
function(qBatch, refType, rem.highCt=FALSE, thresh.Ct=30){
	require(affy)
	
	if( class(refType) == "numeric" ){
		refExp <- qBatch@exprs[,refType]
		refInd <- refType
	}
	else if( refType == "mean" ){
		avg.col <- apply(qBatch@exprs, 2, mean, na.rm=TRUE) 
		refInd <- trunc(median(rank(avg.col)))
		refExp <- qBatch@exprs[,refInd]
	}
	else if( refType == "median" ){
		med.col <- apply(qBatch@exprs, 2, median, na.rm=TRUE) 
		refInd <- trunc(median(rank(avg.col)))
		refExp <- qBatch@exprs[,refInd]
	}
	else if( refType == "pseudo.mean" ){
		refExp <- apply(qBatch@exprs, 1, mean, na.rm=TRUE)
		refInd <- 0 
	}
	else if( refType == "pseudo.median" ){
		refExp <- apply(qBatch@exprs, 1, median, na.rm=TRUE) 
		refInd <- 0 
	}
	
	iv.set <- NULL  
	for( i in 1:ncol(qBatch@exprs) ){
		if( i != refInd ){
			iv.set[[i]] <- as.list(normalize.invariantset(refExp, qBatch@exprs[,i]))
		}
	}
	
	all.genes <- NULL ; numGenes <- length(iv.set[[1]][[2]])
	for( j in 1:length(iv.set) ){
		if( j != refInd ){
			all.genes <- c(all.genes, qBatch@geneNames[iv.set[[j]][[2]]])
		}
	}
	all.genes <- unique(all.genes)
	
	all.genes.ivscore <- integer(length(all.genes)) 
	for( k in 1:length(all.genes) ){
		this.gene <- all.genes[k]
		for( m in 1:length(iv.set) ){
			if( m != refInd ){
				if( this.gene %in% qBatch@geneNames[iv.set[[m]][[2]]] ){
					all.genes.ivscore[k] <- all.genes.ivscore[k] + 1 
				}
			}
		}
	}

	iv.genes <- all.genes[all.genes.ivscore == (length(iv.set)-1)]	
	iv.genes.counts <- NULL 
	
	if( length(iv.genes) == 0 ){
		e <- "No rank invariant genes were found for your data"
		stop(e)
	}
	else{
		for( i in 1:length(iv.genes) ){
			iv.genes.counts[i] <- length(qBatch@geneNames[qBatch@geneNames %in% iv.genes[i]])
		}
		
		avg.iv.expr <- NULL 
		
		for( i in 1:length(iv.genes) ){
			if( iv.genes.counts[i] > 1 ){ 
				xg <- apply(qBatch@exprs[qBatch@geneNames %in% iv.genes[i],], 2, mean)
			}
			else{
				xg <- qBatch@exprs[qBatch@geneNames %in% iv.genes[i],]
			}
			avg.iv.expr <- rbind(avg.iv.expr, xg) 
		}
		
		if( rem.highCt == TRUE ){
			exp.thresh <- NULL 
			for( i in 1:nrow(avg.iv.expr) ){
				exp.thresh <- c(exp.thresh, sum(avg.iv.expr[i,] > thresh.Ct))
			}
			avg.iv.expr <- avg.iv.expr[exp.thresh == 0,]
			iv.genes <- iv.genes[exp.thresh == 0]
		}
		
		# calculate scaling factor 
		
		scale.iv <- apply(avg.iv.expr, 2, mean)
		scale.iv <- scale.iv/scale.iv[1]
			
		snorm.dat <- qBatch@exprs 
		for( i in 1:ncol(snorm.dat) ){
			snorm.dat[,i] <- qBatch@exprs[,i]*scale.iv[i]
		}
			
		colnames(snorm.dat) <- colnames(qBatch@exprs) 
		rownames(snorm.dat) <- qBatch@geneNames 
		qBatch@exprs <- snorm.dat
		qBatch@normalized <- TRUE
		qBatch@normGenes <- as.character(iv.genes)
		return(qBatch) 		
	}
}

`plotVarMean` <-
function(qpcrBatch1, qpcrBatch2, normTag1="Normalization Type1", normTag2="Normalization Type2", ...){
	mainTitle <- paste(normTag1, normTag2, sep=" : ")
	exprs1 <- qpcrBatch1@exprs ; exprs2 <- qpcrBatch2@exprs
	var1 <- apply(exprs1, 1, var, na.rm=TRUE) 
	var2 <- apply(exprs2, 1, var, na.rm=TRUE) 
	var.vals <- log(var1/var2, base=2)
	avg.vals <- apply(cbind(exprs1, exprs2), 1, mean, na.rm=TRUE)
	plot(avg.vals, var.vals, pch=20, xlab="Average Expression", ylab="Log(Variance1:Variance2)", main=mainTitle, ...) 
	abline(h=0, col="blue", lwd=3, lty=2) 
	ind <- is.na(avg.vals) | is.infinite(avg.vals) | is.na(var.vals) | is.infinite(var.vals)
	lines(lowess(avg.vals[!ind], var.vals[!ind]), lwd=3, col="red")
}

`readQpcr` <-
function(fileName, header=TRUE, qc=FALSE, quote="\"", dec=".", fill=TRUE, comment.char="", ...){
	qpcr.data <- read.table(fileName, header=header, quote=quote, dec=dec, fill=fill, comment.char=comment.char, ...)
	primerNames <- as.character(qpcr.data[,1])
	plateID <- as.character(qpcr.data[,2])
	otherNames <- tolower(colnames(qpcr.data)[-(1:2)]) ; whichCt <- grep("ct", otherNames)
	ctVals <- qpcr.data[,whichCt+2] ; ctVals <- data.matrix(ctVals) 
	if( qc ){
		ctVals.qced <- ctQc(ctVals) 
	}
	else{
		ctVals <- apply(ctVals, 1, mean, na.rm=TRUE)
	}
	qres <- new("qpcrBatch", geneNames=primerNames, plateIndex=plateID, exprs = cbind(ctVals), normalized=FALSE, normGenes=primerNames)
	return(qres)
}

`readQpcrBatch` <-
function(..., filenames=character(), qc=FALSE){
	x.data <- NULL
	if(length(filenames) == 0 ){ filenames = dir() }
	for( i in 1:length(filenames) ){
		x.data <- cbind(x.data, readQpcr(filenames[i], qc=qc, ...)@exprs)
		# put something here to check that genes and plates are consistent
	}
	i <- 1 ; x <- readQpcr(filenames[i], qc=qc, ...)
	return(new("qpcrBatch", geneNames = x@primerNames, plateIndex = x@plateIndex, exprs = x.data, normalized = FALSE, normGenes=x@normGenes))
}

`writeQpcr` <-
function(qBatch, fileName, ...){

	xdat <- data.frame(qBatch@geneNames, qBatch@plateIndex, qBatch@exprs)
	colnames(xdat) <- c("GeneNames", "PlateIndex", colnames(qBatch@exprs))
	write.table(xdat, file=fileName, row.names=FALSE, col.names=TRUE)
}




