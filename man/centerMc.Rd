% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxDelaunay.R
\name{centerMc}
\alias{centerMc}
\title{Parameterized center of an interval}
\usage{
centerMc(int, c = 0.5)
}
\arguments{
\item{int}{A \code{vector} with two entries representing an interval.}

\item{c}{A positive real number in \eqn{(0,1)}
parameterizing the center inside \code{int}\eqn{=(a,b)}
with the default \code{c=.5}.
For the interval, \code{int}\eqn{=(a,b)},
the parameterized center is \eqn{M_c=a+c(b-a)}.}
}
\value{
(parameterized) center inside \code{int}
}
\description{
Returns the (parameterized) center, \eqn{M_c},
of the interval, \code{int}\eqn{=(a,b)},
parameterized by \eqn{c \in (0,1)}
so that \eqn{100c} \% of the length of interval is to the left of \eqn{M_c}
and \eqn{100(1-c)} \% of the length of the interval
is to the right of \eqn{M_c}.
That is, for the interval, \code{int}\eqn{=(a,b)},
the parameterized center is \eqn{M_c=a+c(b-a)}.

See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
}
\examples{
c<-.4
a<-0; b<-10
int = c(a,b)
centerMc(int,c)

c<-.3
a<-2; b<-4; int<-c(a,b)
centerMc(int,c)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{centersMc}}
}
\author{
Elvan Ceyhan
}
