\name{qpcrBatch-class}
\docType{class}
\alias{qpcrBatch}
\alias{qpcrBatch-class}
\alias{normalize, qpcrBatch-class}

\title{Class qpcrBatch }
\description{ This is a class representation for qPCR expression data.	}
\section{Objects from the Class}{
Objects can be created using the function \code{\link{readQpcr}} or \code{\link{readQpcrBatch}} to read in 
	raw data from a text file(s). Objects can also be created by using \code{new("qpcrBatch", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{geneNames}:}{Character vector denoting gene or primer pair names. }
    \item{\code{plateIndex}:}{ Character vector denoting plate indices. }
    \item{\code{exprs}:}{ Matrix of qPCR expression values, normally these are the Ct values. }
    \item{\code{normalized}:}{ Logical value, TRUE if expression data has been normalized. }
    \item{\code{normGenes}:}{ Character vector of genes used by the normalization algorithm. }
  }
}
\section{Methods}{
No methods have yet been defined with class "qpcrBatch" in the signature.
}
\author{ Jess Mar \email{jess@jimmy.harvard.edu} }
\note{ This class is better describe in the vignette. }
\examples{
## load example data 
data(qpcrBatch.object)
class(qpcrBatch.object)
}
\keyword{classes}
