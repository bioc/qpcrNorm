\name{matrixByPlate}
\alias{matrixByPlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Internal function to reorganize qPCR data into a rectangular structure. }
\description{
  This function takes a vector of Ct values from a single qPCR experiment and 
  reorganizes it into a matrix structure. Each column in the matrix represents a 
  different plate that was used in the experiment. 
}
\usage{
matrixByPlate(xvec, plateIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xvec}{ Numeric vector of Ct values. }
  \item{plateIndex}{ Character vector, denoting plate index of each gene or primer pair. }
}
\details{
  On the resulting matrix structure: \cr
	The number of rows equals the maximum number of genes or primer pairs that were used on a plate in the experiment. 
	For plates with less genes, NA values are padded at the end of the column vector to complete the rectangular structure. 
	Note: these NA values do not affect downstream calculations. 
}
\value{
  A \code{matrix} object.
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\seealso{ \code{\link{normQpcrQuantile}} }
\examples{
	data(qpcrBatch.object)
	mynormQuant.data <- normQpcrQuantile(qpcrBatch.object) 
}
\keyword{internal}