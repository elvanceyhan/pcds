% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceFunctions.R
\name{ASarc.dens.tri}
\alias{ASarc.dens.tri}
\title{Arc density of Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case}
\usage{
ASarc.dens.tri(Xp, tri, M = "CC", tri.cor = FALSE)
}
\arguments{
\item{Xp}{A set of 2D points which constitute the vertices of the AS-PCD.}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{M}{The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
default is \code{M="CC"} i.e., the circumcenter of \code{tri}.}

\item{tri.cor}{A logical argument for computing the arc density for only the points inside the triangle, \code{tri}
(default is \code{tri.cor=FALSE}), i.e., if \code{tri.cor=TRUE} only the induced digraph with the vertices inside \code{tri} are considered in the
computation of arc density.}
}
\value{
Arc density of AS-PCD whose vertices are the 2D numerical data set, \code{Xp};
AS proximity regions are defined with respect to the triangle \code{tri} and \eqn{CC}-vertex regions.
}
\description{
Returns the arc density of AS-PCD whose vertex set is the given 2D numerical data set, \code{Xp},
(some of its members are) in the triangle \code{tri}.

AS proximity regions is defined with respect to \code{tri}
and vertex regions are defined with the center \code{M="CC"} for circumcenter of \code{tri};
or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
For the number of arcs, loops are not allowed so arcs are only possible for points inside \code{tri}
for this function.

\code{tri.cor} is a logical argument for triangle correction (default is \code{TRUE}),
if \code{TRUE}, only the points inside the triangle are considered (i.e., digraph induced by these vertices
are considered) in computing the arc density, otherwise all points are considered
(for the number of vertices in the denominator of arc density).

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);

set.seed(1)
n<-10  #try also n<-20

Xp<-runif.tri(n,Tr)$g

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)

num.arcsAStri(Xp,Tr,M)
ASarc.dens.tri(Xp,Tr,M)
ASarc.dens.tri(Xp,Tr,M,tri.cor = FALSE)

ASarc.dens.tri(Xp,Tr,M)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{ASarc.dens.tri}}, \code{\link{CSarc.dens.tri}}, and \code{\link{num.arcsAStri}}
}
\author{
Elvan Ceyhan
}
