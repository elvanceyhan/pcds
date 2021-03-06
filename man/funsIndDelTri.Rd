% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{funsIndDelTri}
\alias{funsIndDelTri}
\alias{ind.Del.tri}
\alias{indices.Del.tri}
\title{Functions provide the indices of the Delaunay triangles where the points reside}
\usage{
ind.Del.tri(pt, Yp, DTmesh = NULL)

indices.Del.tri(dat, Yp, DTmesh = NULL)
}
\arguments{
\item{pt}{A 2D point; the index of the Delaunay triangle in which \code{pt} resides is to be
determined. It is an argument for \code{ind.Del.tri}.}

\item{Yp}{A set of 2D points from which Delaunay triangulation is constructed.}

\item{DTmesh}{Delaunay triangles based on \code{Yp}, default is \code{NULL},
which is computed via \code{\link[interp]{tri.mesh}} function
in \code{interp} package. \code{\link[interp]{triangles}} function yields a triangulation data structure from the triangulation object
created by \code{tri.mesh}.}

\item{dat}{A set of 2D points representing the set of data points for which the indices of the Delaunay
triangles they reside is to be determined. It is an argument for \code{indices.Del.tri}.}
}
\value{
ind.Del.tri returns the index of the Delaunay triangle in which the given point, \code{pt}, resides
and indices.Del.tri returns the vector of indices of the Delaunay triangles in which points in the data
set, dat, reside
}
\description{
Two functions: \code{ind.Del.tri} and \code{indices.Del.tri}.

\code{ind.Del.tri} finds the index of the Delaunay triangle in which the given point, \code{pt}, resides.

\code{indices.Del.tri} finds the indices of triangles for all the points in data set, \code{dat}, as a vector.

Delaunay triangulation is based on \code{Yp} and \code{DTmesh} are the Delaunay triangles with default \code{NULL}.
The function returns \code{NA} for a point not inside the convex hull of \code{Yp}.
Number of \code{Yp} points (i.e., size of \code{Yp}) should be at least three and the points
should be in general position so that Delaunay triangulation is (uniquely) defined.

If the number of \code{Yp} points is 3, then there is only one Delaunay triangle and the indices of all
the points inside this triangle are all 1.

See (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay
triangulation and the corresponding algorithm.
}
\examples{
#Examples for ind.Del.tri
nx<-100 #number of X points (target)
ny<-10 #number of Y points (nontarget)

Yp<-cbind(runif(ny),runif(ny))

dat<-runifMT(nx,Yp)$g #data under CSR in the convex hull of Ypoints
#try also dat<-cbind(runif(nx),runif(nx))

ind.Del.tri(dat[10,],Yp)

#or use
DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation
TRY<-interp::triangles(DTY)[,1:3];
ind.Del.tri(dat[10,],Yp,DTY)

ind.Del.tri(c(.5,.5),Yp)

ind.Del.tri(c(1.5,.5),Yp)

ind.DT<-vector()
for (i in 1:nx)
  ind.DT<-c(ind.DT,ind.Del.tri(dat[i,],Yp))
ind.DT

Xlim<-range(Yp[,1],dat[,1])
Ylim<-range(Yp[,2],dat[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation based on Y points

#plot of the data in the convex hull of Y points together with the Delaunay triangulation
plot(dat,main=" ", xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),type="n")
interp::plot.triSht(DTY, add=TRUE, do.points = TRUE,pch=16,col="blue")
points(dat,pch=".",cex=3)
text(dat,labels = factor(ind.DT) )

Yp<-rbind(c(.3,.2),c(.4,.5),c(.14,.15))
ind.Del.tri(c(.25,.25),Yp)

dat.fr<-data.frame(a=Yp)
ind.Del.tri(c(.25,.25),dat.fr)

#Examples for indices.Del.tri
#nx is number of X points (target) and ny is number of Y points (nontarget)
nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;

set.seed(1)
Yp<-cbind(runif(ny),runif(ny))
dat<-runifMT(nx,Yp)$g #data under CSR in the convex hull of Ypoints
#try also dat<-cbind(runif(nx),runif(nx))

tr.ind<-indices.Del.tri(dat,Yp)  #indices of the Delaunay triangles
tr.ind

#or use
DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation based on Y points
tr.ind<-indices.Del.tri(dat,Yp,DTY)  #indices of the Delaunay triangles
tr.ind

Xlim<-range(Yp[,1],dat[,1])
Ylim<-range(Yp[,2],dat[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

#plot of the data in the convex hull of Y points together with the Delaunay triangulation

oldpar <- par(mfrow = c(1,2))
par(pty="s")
plot(dat,main=" ", xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),pch=".")
interp::plot.triSht(DTY, add=TRUE, do.points = TRUE,pch=16,col="blue")
text(dat,labels = factor(tr.ind) )
par(oldpar)

Yp<-rbind(c(.3,.2),c(.4,.5),c(.14,.15))
indices.Del.tri(c(.25,.25),Yp)

dat<-rbind(c(.25,.25),c(.15,.25))
indices.Del.tri(dat,Yp)

dat.fr<-data.frame(a=dat)
indices.Del.tri(dat.fr,Yp)

dat.fr<-data.frame(a=Yp)
indices.Del.tri(c(.25,.25),dat.fr)

}
\references{
\insertAllCited{}
}
