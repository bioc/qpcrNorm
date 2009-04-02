\name{normQpcrHouseKeepingGenes}
\alias{normQpcrHouseKeepingGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for Housekeeping Gene Normalization of qPCR Data. }
\description{
  Implements housekeeping gene normalization for a \code{\link{qpcrBatch}} object.
}
\usage{
normQpcrHouseKeepingGenes(qBatch, hkeep.genes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qBatch}{ A \code{\link{qpcrBatch}} object to be normalized. }
  \item{hkeep.genes}{ Character vector, specifying which housekeeping genes to be used for normalization. }
}
\details{
  The names in \code{hkeep.genes} must be a subset of the gene or primer pair names slot in the \code{\link{qpcrBatch}} object. 
}
\value{
  A \code{\link{qpcrBatch}} object, the \code{normalized} slot is now set at TRUE. 
}
\author{ 
	Jess Mar \email{jess@jimmy.harvard.edu} 
}
\seealso{ 
	\code{\link{normQpcrQuantile}}, \code{\link{normQpcrRankInvariant}}
}
\examples{
	data(qpcrBatch.object)
	mynormHK.data <- normQpcrHouseKeepingGenes(qpcrBatch.object, c("Gpx4"))
}   
\keyword{methods}

