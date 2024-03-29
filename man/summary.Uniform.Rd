% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassFunctions.R
\name{summary.Uniform}
\alias{summary.Uniform}
\title{Return a summary of a \code{Uniform} \code{object}}
\usage{
\method{summary}{Uniform}(object, ...)
}
\arguments{
\item{object}{An \code{object} of class \code{Uniform}.}

\item{\dots}{Additional parameters for \code{summary}.}
}
\value{
The \code{call} of the \code{object} of class \code{"Uniform"},
the \code{type} of the pattern (i.e. the description
of the uniform distribution), study window,
vertices of the support of the Uniform distribution,
some sample points generated from the uniform distribution,
and the number of points (i.e., number of generated
points and the number of vertices of
the support of the uniform distribution.)
}
\description{
Returns the below information about the \code{object}:

\code{call} of the function defining the \code{object},
the \code{type} of the pattern (i.e. the description
of the uniform distribution), study window,
vertices of the support of the Uniform distribution,
some sample points generated from the uniform distribution,
and the number of points (i.e., number of generated
points and the number of vertices of the support
of the uniform distribution.)
}
\examples{
\dontrun{
n<-10 #try also 20, 100, and 1000
A<-c(1,1); B<-c(2,0); R<-c(1.5,2);
Tr<-rbind(A,B,R)

Xdt<-runif.tri(n,Tr)
Xdt
summary(Xdt)
}

}
\seealso{
\code{\link{print.Uniform}},
\code{\link{print.summary.Uniform}},
and \code{\link{plot.Uniform}}
}
