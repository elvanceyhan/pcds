% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentSim2D.R
\name{IarcCSset2pnt.std.tri}
\alias{IarcCSset2pnt.std.tri}
\title{The indicator for the presence of an arc from a point in set \code{S} to the point \code{p} for Central Similarity
Proximity Catch Digraphs (CS-PCDs) - standard equilateral triangle case}
\usage{
IarcCSset2pnt.std.tri(S, p, t, M = c(1, 1, 1))
}
\arguments{
\item{S}{A set of 2D points. Presence of an arc from a point in \code{S} to point \code{p} is checked
by the function.}

\item{p}{A 2D point. Presence of an arc from a point in \code{S} to point \code{p} is checked
by the function.}

\item{t}{A positive real number which serves as the expansion parameter in CS proximity region in the
standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)}
i.e., the center of mass of \eqn{T_e}.}
}
\value{
\eqn{I(}\code{p} is in \eqn{\cup_{x in S} N_{CS}(x,t))}, that is, returns 1 if \code{p} is in \code{S} or inside \eqn{N_{CS}(x,t)} for at least
one \eqn{x} in \code{S}, returns 0 otherwise. CS proximity region is constructed with respect to the standard
equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with \code{M}-edge regions.
}
\description{
Returns \eqn{I(}\code{p} in \eqn{N_{CS}(x,t)} for some \eqn{x} in \code{S}\eqn{)}, that is, returns 1 if \code{p} is in \eqn{\cup_{x in S} N_{CS}(x,t)},
returns 0 otherwise, CS proximity region is constructed with respect to the standard equilateral triangle
\eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with the expansion parameter \eqn{t>0} and edge regions are based
on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e} (which is equivalent to circumcenter of \eqn{T_e}).

Edges of \eqn{T_e}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as edges 3, 1, and 2, respectively.
If \code{p} is not in \code{S} and either \code{p} or all points in \code{S} are outside \eqn{T_e}, it returns 0,
but if \code{p} is in \code{S}, then it always returns 1 regardless of its location (i.e., loops are allowed).

See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
Te<-rbind(A,B,C);
n<-10

set.seed(1)
Xp<-runif.std.tri(n)$gen.points

M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)

t<-.5

S<-rbind(Xp[1,],Xp[2,])  #try also S<-c(.5,.5)
IarcCSset2pnt.std.tri(S,Xp[3,],t,M)
IarcCSset2pnt.std.tri(S,Xp[3,],t=1,M)
IarcCSset2pnt.std.tri(S,Xp[3,],t=1.5,M)

S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
IarcCSset2pnt.std.tri(S,Xp[3,],t,M)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{IarcCSset2pnt.tri}}, \code{\link{IarcCSstd.tri}}, \code{\link{IarcCStri}}, and \code{\link{IarcPEset2pnt.std.tri}}
}
\author{
Elvan Ceyhan
}
