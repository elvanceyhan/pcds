% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassFunctions.R
\name{plot.Lines3D}
\alias{plot.Lines3D}
\title{Plot a \code{Lines3D} \code{object}}
\usage{
\method{plot}{Lines3D}(x, xlab = "x", ylab = "y", zlab = "z", phi = 40, theta = 40, ...)
}
\arguments{
\item{x}{Object of class \code{Lines3D}.}

\item{xlab, ylab, zlab}{Titles for the \eqn{x}, \eqn{y}, and \eqn{z} axes,
respectively (default is \code{xlab="x"}, \code{ylab="y"}
and \code{zlab="z"}).}

\item{theta, phi}{The angles defining the viewing direction.
\code{theta} gives the azimuthal direction and \code{phi} the colatitude.
See \code{\link[plot3D]{persp3D}} for more details.}

\item{\dots}{Additional parameters for \code{plot}.}
}
\value{
None
}
\description{
Plots the line together with the defining vectors
(i.e., the initial and direction vectors).
}
\examples{
\dontrun{
P<-c(1,10,3); Q<-c(1,1,3);
vecs<-rbind(P,Q)
Line3D(P,Q,.1)
Line3D(P,Q,.1,dir.vec=FALSE)

tr<-range(vecs);
tf<-(tr[2]-tr[1])*.1
#how far to go at the lower and upper ends in the x-coordinate
tsq<-seq(-tf*10-tf,tf*10+tf,l=3) #try also l=10, 20 or 100

lnPQ3D<-Line3D(P,Q,tsq)
lnPQ3D
plot(lnPQ3D)
}

}
\seealso{
\code{\link{print.Lines3D}}, \code{\link{summary.Lines3D}},
and \code{\link{print.summary.Lines3D}}
}
