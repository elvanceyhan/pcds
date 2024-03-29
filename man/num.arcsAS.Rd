% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceFunctions.R
\name{num.arcsAS}
\alias{num.arcsAS}
\title{Number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs) and related
quantities of the induced subdigraphs for points in the Delaunay triangles - multiple triangle case}
\usage{
num.arcsAS(Xp, Yp, M = "CC")
}
\arguments{
\item{Xp}{A set of 2D points which constitute the vertices of the AS-PCD.}

\item{Yp}{A set of 2D points which constitute the vertices of the Delaunay triangles.}

\item{M}{The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle or 3D point in barycentric
coordinates which serves as a center in the interior of each Delaunay triangle;
default is \code{M="CC"} i.e., the circumcenter of each triangle.}
}
\value{
A \code{list} with the elements
\item{desc}{A short description of the output: number of arcs
and related quantities for the induced subdigraphs in the Delaunay triangles}
\item{num.arcs}{Total number of arcs in all triangles, i.e., the number of arcs for the entire AS-PCD}
\item{num.in.conhull}{Number of \code{Xp} points in the convex hull of \code{Yp} points}
\item{num.in.tris}{The vector of number of \code{Xp} points in the Delaunay triangles based on \code{Yp} points}
\item{weight.vec}{The \code{vector} of the areas of Delaunay triangles based on \code{Yp} points}
\item{tri.num.arcs}{The \code{vector} of the number of arcs of the component of the AS-PCD in the
Delaunay triangles based on \code{Yp} points}
\item{del.tri.ind}{A matrix of indices of Delaunay triangles based on \code{Yp} points,
each column corresponds to the vector of indices of the vertices of one of the Delaunay triangle.}
\item{data.tri.ind}{A \code{vector} of indices of vertices of the Delaunay triangles in which data points reside,
i.e., column number of \code{del.tri.ind} for each \code{Xp} point.}
\item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
is the Delaunay triangulation based on \code{Yp} points.}
\item{vertices}{Vertices of the digraph, \code{Xp}.}
}
\description{
An object of class \code{"NumArcs"}.
Returns the number of arcs and various other quantities related to the Delaunay triangles
for Arc Slice Proximity Catch Digraph (AS-PCD) whose vertices are the data points in \code{Xp}
in the multiple triangle case.

AS proximity regions are defined with respect to the
Delaunay triangles based on \code{Yp} points and vertex regions in each triangle
are based on the center \code{M="CC"} for circumcenter of each Delaunay triangle
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of each Delaunay triangle;
default is \code{M="CC"} i.e., circumcenter of each triangle.

Convex hull of \code{Yp} is partitioned by the Delaunay triangles based on \code{Yp} points
(i.e., multiple triangles are the set of these Delaunay triangles whose union constitutes
the convex hull of \code{Yp} points).

See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation
and the corresponding algorithm.
}
\examples{
\dontrun{
nx<-15; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;

set.seed(1)
Xp<-cbind(runif(nx),runif(nx))
Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))

M<-"CC"  #try also M<-c(1,1,1)
Narcs = num.arcsAS(Xp,Yp,M)
Narcs
summary(Narcs)
plot(Narcs)

}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{num.arcsAStri}}, \code{\link{num.arcsPE}}, and \code{\link{num.arcsCS}}
}
\author{
Elvan Ceyhan
}
