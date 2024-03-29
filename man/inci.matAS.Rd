% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceFunctions.R
\name{inci.matAS}
\alias{inci.matAS}
\title{Incidence matrix for Arc Slice Proximity Catch Digraphs (AS-PCDs) - multiple triangle case}
\usage{
inci.matAS(Xp, Yp, M = "CC")
}
\arguments{
\item{Xp}{A set of 2D points which constitute the vertices of the AS-PCD.}

\item{Yp}{A set of 2D points which constitute the vertices of the Delaunay triangles.}

\item{M}{The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle or 3D point in barycentric
coordinates which serves as a center in the interior of each Delaunay triangle;
default is \code{M="CC"} i.e., the circumcenter of each triangle.}
}
\value{
Incidence matrix for the AS-PCD whose vertices are the 2D data set, \code{Xp},
and AS proximity regions are defined in the Delaunay triangles based on \code{Yp} points.
}
\description{
Returns the incidence matrix for the AS-PCD whose vertices are a given 2D numerical data set, \code{Xp},
in the convex hull of \code{Yp} which is partitioned by the Delaunay triangles based on \code{Yp} points.

AS proximity regions are defined with respect to the Delaunay triangles based on \code{Yp} points and vertex
regions are based on the center \code{M="CC"}
for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.
Loops are allowed, so the diagonal entries are all equal to 1.

See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation
and the corresponding algorithm.
}
\examples{
\dontrun{
#nx is number of X points (target) and ny is number of Y points (nontarget)
nx<-15; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;

set.seed(1)
Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))

M<-"CC"  #try also M<-c(1,1,1)

IM<-inci.matAS(Xp,Yp,M)
IM
dom.num.greedy(IM)  #try also dom.num.exact(IM)  #this might take a long time for large  nx

IM<-inci.matAS(Xp,Yp[1:3,],M)

inci.matAS(Xp,rbind(Yp,Yp))
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{inci.matAStri}}, \code{\link{inci.matPE}}, and \code{\link{inci.matCS}}
}
\author{
Elvan Ceyhan
}
