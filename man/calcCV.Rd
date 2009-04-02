\name{calcCV}
\alias{calcCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the Average Gene-Specific Coefficient of Variation }
\description{
  This function calculates the coefficient of variation for each gene in the qPCR experiment, 
  and returns the average coefficient of variation across all genes.
}
\usage{
calcCV(qBatch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qBatch}{ A \code{qpcrBatch} object. }
}
\value{
	A \code{numeric} value. 
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\examples{
	data(qpcrBatch.object) 
	mynormRI.data <- normQpcrRankInvariant(qpcrBatch.object, 1)
	mynormQuant.data <- normQpcrQuantile(qpcrBatch.object)    
	barplot(c(calcCV(mynormRI.data), calcCV(mynormQuant.data)), col=c("red", "blue"))
}
\keyword{univar}