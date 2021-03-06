% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{fr2vVRCC}
\alias{fr2vVRCC}
\title{The furthest points in a data set from vertices in each CC-vertex region in a triangle}
\usage{
fr2vVRCC(Dt, tri, ch.all.intri = FALSE)
}
\arguments{
\item{Dt}{A set of 2D points representing the set of data points.}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{ch.all.intri}{A logical argument (default=\code{FALSE}) to check whether all data points are inside
the triangle \code{tri}. So if it is \code{TRUE},
the function checks if all data points are inside the closure of the triangle (i.e. interior and boundary
combined) else it does not.}
}
\value{
A list with the elements
\item{txt1}{Vertex labels are A=1, B=2, and C=3 (corresponds to row number in Extrema Points).}
\item{txt2}{A short description of the distances as "Distances to Vertices...".}
\item{type}{Type of the extrema points}
\item{desc}{A short description of the extrema points}
\item{mtitle}{The "main" title for the plot of the exrema}
\item{ext}{The extrema points, here, furthest points from vertices in each CC-vertex region in the
triangle \code{tri}.}
\item{X}{The input data, \code{Dt}, can be a matrix or data frame}
\item{num.points}{The number of data points, i.e., size of \code{Dt}}
\item{supp}{Support of the data points, here, it is the triangle \code{tri} for this function.}
\item{cent}{The center point used for construction of edge regions.}
\item{ncent}{Name of the center, \code{cent}, it is circumcenter "CC" for this function}
\item{regions}{CC-Vertex regions inside the triangle, \code{tri}, provided as a list}
\item{region.names}{Names of the vertex regions as "vr=1","vr=2","vr=3"}
\item{region.centers}{Centers of mass of the vertex regions inside \code{tri}}
\item{dist2ref}{Distances from furthest points in each vertex region to the corresponding vertex}
}
\description{
An object of class "Extrema".
Returns the furthest data points among the data set, \code{Dt}, in each CC-vertex region from the vertex in the
triangle, \eqn{tri=T(A,B,C)}. Vertex region labels/numbers correspond to the row number of the vertex in \code{tri}.
\code{ch.all.intri} is for checking whether all data points are inside \code{tri} (default is \code{FALSE}).

If some of the data points are not inside \code{tri} and \code{ch.all.intri}=\code{TRUE}, then the function yields
an error message.
If some of the data points are not inside \code{tri} and \code{ch.all.intri}=\code{FALSE}, then the function yields
the closest points to edges among the data points inside \code{tri} (yields \code{NA} if there are no data points
inside \code{tri}).

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10  #try also n<-20

set.seed(1)
dat<-runif.tri(n,Tr)$g

Ext<-fr2vVRCC(dat,Tr)
Ext
summary(Ext)
plot(Ext)

fr2vVRCC(dat[1,],Tr)
f2v<-fr2vVRCC(dat,Tr)

CC<-circ.cent.tri(Tr)  #the circumcenter
D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
Ds<-rbind(D1,D2,D3)

Xlim<-range(Tr[,1],dat[,1])
Ylim<-range(Tr[,2],dat[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(Tr,xlab="",asp=1,ylab="",pch=".",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty=2)
points(dat)
points(rbind(f2v$ext),pch=4,col=2)

txt<-rbind(Tr,CC,Ds)
xc<-txt[,1]+c(-.06,.08,.05,.12,-.1,-.1,-.09)
yc<-txt[,2]+c(.02,-.02,.05,.0,.02,.06,-.04)
txt.str<-c("A","B","C","CC","D1","D2","D3")
text(xc,yc,txt.str)

fr2vVRCC(dat,Tr)

fr2vVRCC(c(1.4,1.2),Tr)

dat.fr<-data.frame(a=dat)
fr2vVRCC(dat.fr,Tr)

dat.fr<-data.frame(a=Tr)
fr2vVRCC(dat,dat.fr)

dat2<-rbind(dat,c(.2,.4))
fr2vVRCC(dat2,Tr)
\dontrun{
fr2vVRCC(dat2,Tr,ch.all.intri = TRUE)
#gives an error message since not all points in the data set are in the triangle
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{fr2vTbVRCC}}, \code{\link{fr2eTeER}}, \code{\link{Kfr2vTbVRCC}} and \code{\link{Kfr2vVRCC}}
}
