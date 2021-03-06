% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{NumArcsAStri}
\alias{NumArcsAStri}
\title{Number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case}
\usage{
NumArcsAStri(dat, tri, M = "CC")
}
\arguments{
\item{dat}{A set of 2D points which constitute the vertices of the digraph (i.e., AS-PCD).}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{M}{The center of the triangle. "CC" stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
default is \code{M}="CC" i.e. the circumcenter of \code{tri}.}
}
\value{
The number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs)
}
\description{
Returns the number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs)
whose vertices are the 2D data set, \code{dat}. The data points could be inside or outside a general
triangle \eqn{tri=T(A,B,C)=(rv=1,rv=2,rv=3)}, with vertices of \code{tri} stacked row-wise.

AS proximity regions are defined with respect to the triangle \code{tri} and vertex regions are
based on the center \code{M}="CC" for circumcenter of \code{tri};
or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of the triangle \code{tri}; default is \code{M}="CC" i.e. circumcenter of \code{tri}.
For the number of arcs, loops are not allowed,
so arcs are only possible for points inside the triangle, \code{tri}.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);

n<-10  #try also n<-20
set.seed(1)
dat<-runif.tri(n,Tr)$g

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)

NumArcsAStri(dat,Tr)  #with default M="CC"
NumArcsAStri(dat,Tr,M)

NumArcsAStri(rbind(dat,c(0,2)),Tr,M)

dat.fr<-data.frame(a=dat)
NumArcsAStri(dat.fr,Tr,M)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{NumArcsASMT}}, \code{\link{NumArcsPEtri}}, and \code{\link{NumArcsCStri}}
}
