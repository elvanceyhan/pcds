% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PatternGen.R
\name{runif.std.tri}
\alias{runif.std.tri}
\title{Generation of Uniform Points in the Standard Equilateral Triangle}
\usage{
runif.std.tri(n)
}
\arguments{
\item{n}{A positive integer representing the number of uniform points
to be generated in the standard equilateral triangle \eqn{T_e}.}
}
\value{
A \code{list} with the elements
\item{type}{The type of the pattern from which points are to be generated}
\item{mtitle}{The \code{"main"} title for the plot of the point pattern}
\item{tess.points}{The vertices of the support region of
the uniformly generated points, it is the
standard equilateral triangle \eqn{T_e} for this function}
\item{gen.points}{The output set of generated points uniformly
in the standard equilateral triangle \eqn{T_e}.}
\item{out.region}{The outer region which contains the support region,
\code{NULL} for this function.}
\item{desc.pat}{Description of the point pattern
from which points are to be generated}
\item{num.points}{The \code{vector} of two numbers,
which are the number of generated points and the number
of vertices of the support points (here it is 3).}
\item{txt4pnts}{Description of the two numbers in \code{num.points}}
\item{xlimit,ylimit}{The ranges of the \eqn{x}-
and \eqn{y}-coordinates of the support, \eqn{T_e}}
}
\description{
An object of class \code{"Uniform"}.
Generates \code{n} points uniformly
in the standard equilateral triangle \eqn{T_e=T(A,B,C)}
with vertices \eqn{A=(0,0)}, \eqn{B=(1,0)}, and \eqn{C=(1/2,\sqrt{3}/2)}.
}
\examples{
\dontrun{
A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
Te<-rbind(A,B,C);
n<-100

set.seed(1)
Xdt<-runif.std.tri(n)
Xdt
summary(Xdt)
plot(Xdt,asp=1)

Xlim<-range(Te[,1])
Ylim<-range(Te[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

Xp<-runif.std.tri(n)$gen.points
plot(Te,asp=1,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.01,.01),
ylim=Ylim+yd*c(-.01,.01))
polygon(Te)
points(Xp)
}

}
\seealso{
\code{\link{runif.basic.tri}}, \code{\link{runif.tri}},
and \code{\link{runif.multi.tri}}
}
\author{
Elvan Ceyhan
}
