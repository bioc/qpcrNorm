\name{plotVarMean}
\alias{plotVarMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Constructs scatter plot to compare the effects of two normalization algorithms on a qPCR dataset. }
\description{
  This function makes a scatter plot which serves as a useful exploratory tool in evaluating whether one 
  normalization algorithm has been more effective than another on a given qPCR dataset.
}
\usage{
plotVarMean(qpcrBatch1, qpcrBatch2, normTag1 = "Normalization Type1", normTag2 = "Normalization Type2", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qpcrBatch1}{ A \code{\link{qpcrBatch}} object. }
  \item{qpcrBatch2}{ A \code{\link{qpcrBatch}} object. }
  \item{normTag1}{ Character string denoting what normalization algorithm was used for this data set. }
  \item{normTag2}{ Character string denoting what normalization algorithm was used for this data set. }
  \item{\dots}{ Further arguments can be supplied to the \code{\link{plot}} function. }
}
\details{
  For each gene, the function plots its log-transformed ratio of its expression variance in 
  one normalized dataset versus another normalized dataset, i.e. let Gij be the variance of 
  the expression values of gene i that have been normalized with method j. 
  We plot the natural log-transformed ratio of Gij to Gik on the y-axis, 
  and the average expression of gene i on the x-axis for all genes. /cr
  The red curve represents a smoothed lowess curve that has been fitted to reflect the overall 
  trend of the data. When the red curve drops below y = 0 (the blue dotted line) we know that method j 
  effects a greater reduction in the variation of the data over method k. 
  Similarly, when the red curve is above y = 0, method k is more effective in reducing the variation 
  in the data than method j. If the data from both methods have similar variances then the red curve 
  should remain at y = 0. Bolstad et al. (2003) originally used these plots for variance comparisons of 
  different normalization methods for high density oligonucleotide array data.
}
\value{
  A \code{plot} object. 
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\seealso{ \code{\link{plot}} }
\references{
	Bolstad B et al. A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics, 2003. 
	}
\examples{
  # data(qpcrBatch.object)
  # mynormRI.data <- normQpcrRankInvariant(qpcrBatch.object, 1) 
  # mynormQuant.data <- normQpcrQuantile(qpcrBatch.object)
  # plotVarMean(mynormRI.data, mynormQuant.data, normTag1="Rank-Invariant", normTag2="Quantile", main="Comparing Two Data-driven Methods")  
}
\keyword{aplot}