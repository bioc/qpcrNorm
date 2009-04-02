\name{ctQc}
\alias{ctQc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quality Control Filter for Replicate Ct Values }
\description{
  This function applies a quality control filter to triplicate Ct values before 
  combining them into a single summary Ct measure. \cr
  
  Current implementation can only handle three replicates. 
}
\usage{
ctQc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Matrix with three columns, corresponding to the triplicate Ct values. }
}
\details{
  Applying ctQc is an alternative to averaging the triplicate Ct values. This filter was originally developed by Yasumasa Kimura 
  1. For each primer, sort Ct values in ascending order so we have [ct1, ct2, ct3].
  2. Take two differences 
	difference1 := ct2 - ct1
	difference2 := ct3 - ct2
  If either or both have ct 40, we don't calculate the difference and set it to region 4 in the next step.
  3. Take an average according to the differences.
  Here we have 2 thresholds : 0.2 and 1.0. With the thresholds, we classify the 2 differences into the below regions.
    region1 : difference <= 0.2
    region2 : 0.2 < difference <= 1.0
    region3 : 1.0 < difference
    region4 : either or both ct are 40
  If the 2 differences are in same region, we take an average of 3 ct values.
  If the 2 differences are in different regions, we take an average of 2 ct values which are in smaller number region.
 }
\value{
  Numeric vector of Ct values combined over the three replicates according to the QC filter.
}
\author{ Yasumasa Kimura \cr
		 Jess Mar \email{jess@jimmy.harvard.edu} 
}
\seealso{
	\code{\link{readQpcr}}, \code{\link{readQpcrBatch}}
}
\examples{
	## myQpcrBatch <- readQpcrBatch(qc=T)	# reads in data from a batch of qPCR experiments, applies ctQc to raw Ct values
}
\keyword{internal}