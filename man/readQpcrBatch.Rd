\name{readQpcrBatch}
\alias{readQpcrBatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data Input Function for a Batch of qPCR Experiments. }
\description{
  This function reads in data from multiple qPCR experiments from the one batch. 
  Each text file in the batch must meet the structure required by \code{\link{readQpcr}}. \cr
  Note: In order to qualify as a batch, it is assumed that the same set of primers 
  are being analyzed in each experiment.
}
\usage{
readQpcrBatch(..., filenames = character(), qc = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{ Filenames separated by a comma. }
  \item{filenames}{ Character vector specifying file names. }
  \item{qc}{ Logical value, TRUE if a QC filter \code{\link{ctQc}} should be applied to the data. \cr
				If qc = F, the replicate Ct values will be averaged. See \code{\link{ctQc}}. }
}
\details{
If the function is called with no arguments \code{readQpcrBatch()} all the files in the working directory are 
read and put into a \code{\link{qpcrBatch}} object. 
All files must conform to the following structure: \cr
 	1st column = names denoting genes or primer pairs \cr
	2nd column = plate index of each gene or primer pair \cr
	remaining columns = (replicate) Ct values \cr
	other annotation fields may be added after the Ct value columns. /cr /cr

Currently fields appearing after the Ct value columns are not associated with the output of this function.
  Note: the majority of arguments to readQpcr are identical to those supplied to read.table. These have been included to 
  give the user greater control over data input, should the data deviate from a standard tab-delimited file structure. 
  For a set of standard tab-delimited text files, specifying the \code{filenames} should be sufficient.   
}
\value{
	A \code{\link{qpcrBatch}} object.
}
\examples{
	## myBatch <- readQpcrBatch()
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\seealso{ \code{\link{ctQc}}, \code{\link{readQpcr}}, \code{\link{setwd}} }
\keyword{IO}