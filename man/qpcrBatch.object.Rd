\name{qpcrBatch.object}
\alias{qpcrBatch.object}
\docType{data}
\title{ qpcrBatch instance qpcrBatch.object }
\description{
  This is an artifically generated qPCR data set. The data set 
  has been closely simulated from original data for 2396 genes on 13 time points. 
  Each measurement within the one sample was repeated over three replicate wells, 
  across multiple plates.  
}
\usage{data(qpcrBatch.object)}
\format{
  A data frame with 2396 observations on the following 41 variables.
  \describe{
    \item{\code{Primers}}{ Character vector of gene or primer pair names. } 
    \item{\code{Plate_Index}}{ Numeric vector denoting plate indices. }
    \item{\code{Time1_Rep1}}{ Ct values for first time point, first replicate. }
    \item{\code{Time1_Rep2}}{ Ct values for first time point, second replicate. }
    \item{\code{Time1_Rep3}}{ Ct values for first time point, third replicate. }
  }
}
\examples{
data(qpcrBatch.object)
}
\keyword{datasets}
