% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceFunctions.R
\name{plotASregs}
\alias{plotASregs}
\title{The plot of the Arc Slice (AS) Proximity Regions for a 2D data set - multiple triangle case}
\usage{
plotASregs(
  Xp,
  Yp,
  M = "CC",
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  ...
)
}
\arguments{
\item{Xp}{A set of 2D points for which AS proximity regions are constructed.}

\item{Yp}{A set of 2D points which constitute the vertices of the Delaunay triangulation. The Delaunay
triangles partition the convex hull of \code{Yp} points.}

\item{M}{The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle or 3D point in barycentric
coordinates which serves as a center in the interior of each Delaunay triangle;
default is \code{M="CC"} i.e., the circumcenter of each triangle.}

\item{main}{An overall title for the plot (default=\code{NULL}).}

\item{xlab, ylab}{Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).}

\item{xlim, ylim}{Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
(default=\code{NULL} for both).}

\item{\dots}{Additional \code{plot} parameters.}
}
\value{
Plot of the \code{Xp} points, Delaunay triangles based on \code{Yp} and also the AS proximity regions
for \code{Xp} points inside the convex hull of \code{Yp} points
}
\description{
Plots the \code{Xp} points in and outside of the convex hull of \code{Yp} points and also plots the AS proximity regions
for \code{Xp} points and Delaunay triangles based on \code{Yp} points.

AS proximity regions are constructed with respect
to the Delaunay triangles based on \code{Yp} points (these triangles partition the convex hull of \code{Yp} points),
i.e., AS proximity regions are only defined for \code{Xp} points inside the convex hull of \code{Yp} points.

Vertex regions are based on the center \code{M="CC"}
for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.

See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
}
\examples{
\dontrun{
nx<-10 ; ny<-5

set.seed(1)
Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))

M<-c(1,1,1)  #try also M<-c(1,2,3) #or M="CC"

plotASregs(Xp,Yp,M,xlab="",ylab="")

plotASregs(Xp,Yp[1:3,],M,xlab="",ylab="")

Xp<-c(.5,.5)
plotASregs(Xp,Yp,M,xlab="",ylab="")
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{plotASregs.tri}}, \code{\link{plotPEregs.tri}}, \code{\link{plotPEregs}},
\code{\link{plotCSregs.tri}}, and \code{\link{plotCSregs}}
}
\author{
Elvan Ceyhan
}
