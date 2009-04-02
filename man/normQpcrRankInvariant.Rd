\name{normQpcrRankInvariant}
\alias{normQpcrRankInvariant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for Rank-Invariant Set Normalization for qPCR Data. }
\description{
 Implements rank-invariant set normalization for a \code{\link{qpcrBatch}} object. 
 We have adapted this algorithm from the function \code{\link[affy]{normalize.invariantset}}  
 from the \pkg{affy} package. 
}
\usage{
normQpcrRankInvariant(qBatch, refType, rem.highCt = FALSE, thresh.Ct = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qBatch}{ A \code{\link{qpcrBatch}} object. }
  \item{refType}{ Indicates what reference sample should be used, can be an integer or character string. See Details below. }
  \item{rem.highCt}{ Logical indicator, TRUE if user wishes to remove genes with high Ct values (very low expression) 
						that may be associated poor data quality. }
  \item{thresh.Ct}{ Numerical value indicating the Ct value cutoff threshold, if \code{rem.highCt} = FALSE, 
					genes with Ct values > thresh.Ct are removed from the data set.   }
}
\details{
  The algorithm computes all rank-invariant sets of genes between pairwise comparisons where each 
  experimental sample in the \code{qpcrBatch} object is paired against a reference. There are several ways to specify 
  what a sensible choice for the reference sample should be. \cr
	\cr
  1. The reference is an experimental sample in the \code{qpcrBatch} object. \cr
	 Specify \code{refType} as an \code{integer} value, corresponding to the index of which experimental sample is the reference.  \cr
	\cr
  2. The reference is the sample which is closest to mean of all the experiments. \cr
     Specify \code{refType = "mean"}. \cr
	\cr
  3. The reference is the sample which is closest to median of all the experiments. \cr
     Specify \code{refType = "median"}. \cr
    \cr
  4. The reference is the mean of all experiments in the \code{qpcrBatch} object. \cr
     Specify \code{refType = "pseudo.mean"}. \cr 
    \cr
  5. The reference is the median of all experiments in the \code{qpcrBatch} object. \cr
     Specify \code{refType = "pseudo.median"}. \cr
}
\value{
  A \code{\link{qpcrBatch}} object, the \code{normalized} slot is now set at TRUE. 
  The names of the rank-invariant genes used for normalization are stored as a vector in the \code{normGenes} slot of the \code{qpcrBatch} object returned.
  To retrieve the rank-invariant gene names, use \code{qpcrBatch@normGenes}.
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\seealso{ 
	\code{\link{normQpcrQuantile}}, \code{\link[affy]{normalize.invariantset}}  
}
\examples{
 data(qpcrBatch.object)
 mynormRI.data <- normQpcrRankInvariant(qpcrBatch.object, 1) 
 mynormRI.data@normGenes		# retrieves names of genes in the rank-invariant set
}
\keyword{methods}