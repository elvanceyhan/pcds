% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PropEdge2D.R
\name{Idom.setPEtri}
\alias{Idom.setPEtri}
\title{The indicator for the set of points \code{S} being a dominating set
or not for Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
one triangle case}
\usage{
Idom.setPEtri(S, Xp, tri, r, M = c(1, 1, 1))
}
\arguments{
\item{S}{A set of 2D points which is to be tested for
being a dominating set for the PE-PCDs.}

\item{Xp}{A set of 2D points
which constitute the vertices of the PE-PCD.}

\item{tri}{A \eqn{3 \times 2} matrix with each row
representing a vertex of the triangle.}

\item{r}{A positive real number
which serves as the expansion parameter in PE proximity region
constructed in the triangle \code{tri}; must be \eqn{\ge 1}.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri}
or the circumcenter of \code{tri}
which may be entered as "CC" as well;
default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.}
}
\value{
\eqn{I(}\code{S} a dominating set of PE-PCD\eqn{)},
that is, returns 1 if \code{S} is a dominating set of PE-PCD whose
vertices are the data points in \code{Xp}; and returns 0 otherwise,
where PE proximity region is constructed in
the triangle \code{tri}.
}
\description{
Returns \eqn{I(}\code{S} a dominating set of PE-PCD
whose vertices are the data set \code{Xp}\eqn{)}, that is,
returns 1 if \code{S} is a dominating set of PE-PCD,
and returns 0 otherwise.

PE proximity region is constructed with
respect to the triangle \code{tri}
with the expansion parameter \eqn{r \ge 1} and vertex regions are based
on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
in the interior of the triangle \code{tri}
or based on the circumcenter of \code{tri};
default is \eqn{M=(1,1,1)}, i.e.,
the center of mass of \code{tri}.
The triangle \code{tri}\eqn{=T(A,B,C)} has edges \eqn{AB}, \eqn{BC}, \eqn{AC}
which are also labeled as edges 3, 1, and 2, respectively.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10

set.seed(1)
Xp<-runif.tri(n,Tr)$gen.points

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)

r<-1.5

S<-rbind(Xp[1,],Xp[2,])
Idom.setPEtri(S,Xp,Tr,r,M)

S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
Idom.setPEtri(S,Xp,Tr,r,M)

S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
Idom.setPEtri(S,Xp,Tr,r,M)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{Idom.setPEstd.tri}}, \code{\link{IarcPEset2pnt.tri}},
\code{\link{Idom.setCStri}}, and \code{\link{Idom.setAStri}}
}
\author{
Elvan Ceyhan
}
