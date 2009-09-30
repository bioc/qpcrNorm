\name{normQpcrQuantile}
\alias{normQpcrQuantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for Quantile Normalization of qPCR Data. }
\description{
	Implements quantile normalization for a \code{\link{qpcrBatch}} object. 
	We have adapted this algorithm from the function
	\code{\link[limma:normalizebetweenarrays]{normalizeBetweenArrays}}
	from the \pkg{limma} package.  
	Data in a \code{\link{qpcrBatch}} object is normalized such that
	within an experiment, the expression distributions  
	across plates are more or less identical, and across
	experiments, the expression distributions  
	are also now more or less identical. 
}
\usage{
normQpcrQuantile(qBatch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qBatch}{ A \code{link{qpcrBatch}} object. }
}
\value{
  A \code{link{qpcrBatch}} object, the \code{normalized} slot is now set at TRUE.
}
\author{ 
	Jess Mar \email{jess@jimmy.harvard.edu} 
}
\seealso{ 
  \code{\link{normQpcrRankInvariant}},
  \code{\link[limma:normalizebetweenarrays]{normalizeBetweenArrays}} 
}
\examples{
 data(qpcrBatch.object) 
 mynormQuant.data <- normQpcrQuantile(qpcrBatch.object) 
}
\keyword{methods}
