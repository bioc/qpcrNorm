\name{readQpcr}
\alias{readQpcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data Input Function for a Single qPCR Experiment. }
\description{
  This function reads in data from a single qPCR experiment. The text file must have the following structure: \cr
	1st column = names denoting genes or primer pairs \cr
	2nd column = plate index of each gene or primer pair \cr
	remaining columns = (replicate) Ct values. \cr
	Note: other annotation fields may be added after the Ct value columns. \cr
}
\usage{
readQpcr(fileName, header = TRUE, qc = FALSE, quote = "\"", dec = ".", fill = TRUE, comment.char = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{ Character string. }
  \item{header}{ Logical value, TRUE if the file contains the names of the variables as its first line. }
  \item{qc}{ Logical value, TRUE if a QC filter \code{\link{ctQc}} should be applied to the data. \cr
				If qc = F, the replicate Ct values will be averaged. }
  \item{quote}{ Set of quoting characters. To disable quoting, set quote = "". See \code{\link{scan}} for behaviour on quotes embedded in quotes. }
  \item{dec}{ Character used for decimal points. }
  \item{fill}{ Logical value, TRUE if in case rows have unequal length, blank fields are implicitly added. See \code{\link{read.table}}. }
  \item{comment.char}{ Character vector of length one containing a single character or an empty string. Use "" to turn off the interpretation of comments altogether.  }
  \item{\dots}{ further arguments to be passed to \code{\link{read.table}}. }
}
\details{
  Currently fields appearing after the Ct value columns are not associated with the output of this function.
  Note: the majority of arguments to readQpcr are identical to those supplied to read.table. These have been included to 
  give the user greater control over data input, should the data deviate from a standard tab-delimited file structure. 
  For a standard tab-delimited text file, specifying the \code{fileName} should be sufficient.   
}
\value{
  A \code{\link{qpcrBatch}} object. 
} 

\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\seealso{ \code{\link{readQpcrBatch}}, \code{\link{ctQc}} }

\examples{ 
	## onerun.data <- readQpcr("singleQpcrRun.txt")
}
\keyword{IO}


boo sarita rocks