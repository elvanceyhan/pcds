% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{IndCSTeSet}
\alias{IndCSTeSet}
\title{The indicator for the presence of an arc from a point in set \code{S} to the point \code{pt} for Central Similarity
Proximity Catch Digraphs (CS-PCDs) - standard equilateral triangle case}
\usage{
IndCSTeSet(S, pt, t, M = c(1, 1, 1))
}
\arguments{
\item{S}{A set of 2D points. Presence of an arc from a point in \code{S} to point \code{pt} is checked
by the function.}

\item{pt}{A 2D point. Presence of an arc from a point in \code{S} to point \code{pt} is checked
by the function.}

\item{t}{A positive real number which serves as the expansion parameter in CS proximity region in the
standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)}
i.e. the center of mass of \eqn{T_e}.}
}
\value{
I(\code{pt} is in \eqn{\cup_{x in S} NCS(x,t)}), that is, returns 1 if \code{pt} is in \code{S} or inside \eqn{NCS(x,t)} for at least
one x in \code{S}, returns 0 otherwise. CS proximity region is constructed with respect to the standard
equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with \code{M}-edge regions.
}
\description{
Returns I(\code{pt} in \eqn{NCS(x,t)} for some x in \code{S}), that is, returns 1 if \code{pt} is in \eqn{\cup_{x in S} NCS(x,t)},
returns 0 otherwise, CS proximity region is constructed with respect to the standard equilateral triangle
\eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with the expansion parameter \eqn{t>0} and edge regions are based
on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e. the center of mass of \eqn{T_e} (which is equivalent to circumcenter of \eqn{T_e}).

Edges of \eqn{T_e}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as edges 3, 1, and 2, respectively.
If \code{pt} is not in \code{S} and either \code{pt} or all points in \code{S} are outside \eqn{T_e}, it returns 0,
but if \code{pt} is in \code{S}, then it always returns 1 regardless of its location (i.e., loops are allowed).

See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
}
\examples{
A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
Te<-rbind(A,B,C);
n<-10

set.seed(1)
dat<-runifTe(n)$gen.points

M<-as.numeric(runifTe(1)$g)  #try also M<-c(.6,.2)

t<-.5

S<-rbind(dat[1,],dat[2,])  #try also S<-c(.5,.5)
IndCSTeSet(S,dat[3,],t,M)
IndCSTeSet(S,dat[3,],t=1,M)
IndCSTeSet(S,dat[3,],t=1.5,M)

S<-rbind(dat[1,],dat[2,],dat[3,],dat[5,])
IndCSTeSet(S,dat[3,],t,M)

IndCSTeSet(S,dat[6,],t,M)
IndCSTeSet(S,dat[6,],t=.25,M)

S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
IndCSTeSet(S,dat[3,],t,M)

IndCSTeSet(c(.2,.5),dat[2,],t,M)
IndCSTeSet(dat,c(.2,.5),t,M)
IndCSTeSet(dat,dat[2,],t)
IndCSTeSet(c(.2,.5),c(.2,.5),t,M)
IndCSTeSet(dat[5,],dat[2,],t,M)

S<-rbind(dat[1,],dat[2,],dat[3,],dat[5,],c(.2,.5))
IndCSTeSet(S,dat[3,],t,M)

P<-c(.4,.2)
S<-dat[c(1,3,4),]
IndCSTeSet(dat,P,t,M)

IndCSTeSet(rbind(S,S),P,t,M)

dat.fr<-data.frame(a=S)
IndCSTeSet(dat.fr,P,t,M)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{IndNCStriSet}}, \code{\link{IndCSTe}}, \code{\link{IndNCStri}}, and \code{\link{IndNPETeSet}}
}
