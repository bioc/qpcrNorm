\name{writeQpcr}
\alias{writeQpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Writes qpcrBatch object out to a File. }
\description{
  This function writes a \code{\link{qpcrBatch}} out to a tab-delimited text file. 
  \code{\link{writeQpcr}} can be used to write out the normalized qPCR data out 
  to an external file.  
}
\usage{
  writeQpcr(qBatch, fileName, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qBatch}{ A \code{\link{qpcrBatch}} object. }
  \item{fileName}{ Character string specifying name of the output file. }
  \item{\dots}{ Extra arguments to be passed to \code{\link{write.table}}. }
}
\details{
  Function creates a tab-delimited text file with three columns, \cr 
	1st column = names denoting genes or primer pairs 
	2nd column = plate index
	3rd column = normalized Ct value
}
\references{ 
Mar J et al. Data-driven Normalization Strategies for qPCR Data. Technical Report, 2008.
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\seealso{ \code{\link{write.table}} }
\examples{
	## writeQpcr(qpcrBatch.object, "output1.txt")
}
\keyword{IO}