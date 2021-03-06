% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{IndNAStriSet}
\alias{IndNAStriSet}
\title{The indicator for the presence of an arc from a point in set \code{S} to the point \code{pt} for
Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case}
\usage{
IndNAStriSet(S, pt, tri, M = "CC")
}
\arguments{
\item{S}{A set of 2D points whose AS proximity regions are considered.}

\item{pt}{A 2D point. The function determines whether \code{pt} is inside the union of AS proximity
regions of points in \code{S} or not.}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{M}{The center of the triangle. "CC" stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
default is \code{M}="CC" i.e. the circumcenter of \code{tri}.}
}
\value{
I(\code{pt} in U_{x in \code{S}}NAS(x,r)), that is, returns 1 if \code{pt} is in \code{S} or inside NAS(x) for at least
one x in \code{S}, returns 0 otherwise, where AS proximity region is constructed in tri
}
\description{
Returns I(\eqn{pt \in N_{AS}(x)} for some \eqn{x \in S}), that is, returns 1 if \code{pt} is in \eqn{\cup_{x \in S}NAS(x)},
returns 0 otherwise, where \eqn{N_{AS}(x)} is the AS proximity region for point \eqn{x}.

AS proximity regions are constructed with respect to the triangle, \eqn{tri=T(A,B,C)=(rv=1,rv=2,rv=3)},
and vertices of \code{tri} are also labeled as 1,2, and 3, respectively.

Vertex regions are based on the center \code{M}="CC" for circumcenter of \code{tri};
or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of the triangle \code{tri}; default is \code{M}="CC" i.e. circumcenter of \code{tri}.

If \code{pt} is not in \code{S} and either \code{pt} or all points in \code{S} are outside \code{tri}, it returns 0,
but if \code{pt} is in \code{S}, then it always returns 1 (i.e., loops are allowed).

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10

set.seed(1)
dat<-runif.tri(n,Tr)$gen.points

S<-rbind(dat[1,],dat[2,])  #try also S<-c(1.5,1)

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)

IndNAStriSet(S,dat[3,],Tr,M)

S<-rbind(dat[1,],dat[2,],dat[3,],dat[5,])
IndNAStriSet(S,dat[3,],Tr,M)

IndNAStriSet(S,dat[6,],Tr,M)

S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
IndNAStriSet(S,dat[3,],Tr,M)

IndNAStriSet(c(.2,.5),dat[2,],Tr,M)
IndNAStriSet(dat,c(.2,.5),Tr,M)
IndNAStriSet(dat,dat[2,],Tr,M)
IndNAStriSet(c(.2,.5),c(.2,.5),Tr,M)
IndNAStriSet(dat[5,],dat[2,],Tr,M)

S<-rbind(dat[1,],dat[2,],dat[3,],dat[5,],c(.2,.5))
IndNAStriSet(S,dat[3,],Tr,M)

P<-c(.4,.2)
S<-dat[c(1,3,4),]
IndNAStriSet(dat,P,Tr,M)
IndNAStriSet(S,P,Tr,M)

IndNAStriSet(rbind(S,S),P,Tr,M)

dat.fr<-data.frame(a=S)
IndNAStriSet(dat.fr,P,Tr,M)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{IndNAStri}}, \code{\link{IndNAStriSet}}, and \code{\link{IndNCStriSet}}
}
