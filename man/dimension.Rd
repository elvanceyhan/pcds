% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxGeometry.R
\name{dimension}
\alias{dimension}
\title{The dimension of a \code{vector} or matrix or a data frame}
\usage{
dimension(x)
}
\arguments{
\item{x}{A \code{vector} or a matrix or a data frame
whose dimension is to be determined.}
}
\value{
Dimension (i.e., number of columns) of \code{x}
}
\description{
Returns the dimension (i.e., number of columns) of \code{x},
which is a matrix or a \code{vector} or a data frame.
This is different than the \code{dim} function in base \code{R},
in the sense that,
\code{dimension} gives only the number of columns of the argument \code{x},
while \code{dim} gives the number of rows and
columns of \code{x}.
\code{dimension} also works for a scalar or a vector,
 while \code{dim} yields \code{NULL} for such arguments.
}
\examples{
\dontrun{
dimension(3)
dim(3)

A<-c(1,2)
dimension(A)
dim(A)

B<-c(2,3)
dimension(rbind(A,B,A))
dimension(cbind(A,B,A))

M<-matrix(runif(20),ncol=5)
dimension(M)
dim(M)

dimension(c("a","b"))
}

}
\seealso{
\code{\link{is.point}} and \code{\link{dim}}
from the base distribution of \code{R}
}
\author{
Elvan Ceyhan
}
