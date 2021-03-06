% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{Kfr2vVRCC}
\alias{Kfr2vVRCC}
\title{The \code{k} furthest points in a data set from vertices in each CC-vertex region in a triangle}
\usage{
Kfr2vVRCC(Dt, tri, k, ch.all.intri = FALSE)
}
\arguments{
\item{Dt}{A set of 2D points representing the set of data points.}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{k}{A positive integer. \code{k} furthest data points in each CC-vertex region are to be found if exists, else
\code{NA} are provided for (some of) the \code{k} furthest points.}

\item{ch.all.intri}{A logical argument (default=\code{FALSE}) to check whether all data points are inside
the triangle \code{tri}. So if it is \code{TRUE},
the function checks if all data points are inside the closure of the triangle (i.e. interior and boundary
combined) else it does not.}
}
\value{
A list with the elements
\item{txt1}{Vertex labels are A=1, B=2, and C=3 (corresponds to row number in Extrema Points).}
\item{txt2}{A shorter description of the distances as "Distances of \code{k} furthest points in the vertex regions
 to Vertices".}
\item{type}{Type of the extrema points}
\item{desc}{A short description of the extrema points}
\item{mtitle}{The "main" title for the plot of the exrema}
\item{ext}{The extrema points, here, \code{k} furthest points from vertices in each CC-vertex region in
the triangle \code{tri}.}
\item{X}{The input data, \code{Dt}, can be a matrix or data frame}
\item{num.points}{The number of data points, i.e., size of \code{Dt}}
\item{supp}{Support of the data points, it is \code{tri} for this function.}
\item{cent}{The center point used for construction of vertex regions}
\item{ncent}{Name of the center, \code{cent}, it is circumcenter "CC" for this function.}
\item{regions}{Vertex regions inside the triangle, \code{tri}, provided as a list}
\item{region.names}{Names of the vertex regions as "vr=1","vr=2","vr=3"}
\item{region.centers}{Centers of mass of the vertex regions inside \eqn{T_b}.}
\item{dist2ref}{Distances from \code{k} furthest points in each vertex region to the corresponding vertex
(each row representing a vertex in \code{tri}). Among the distances the first \code{k} entries are the distances
from the \code{k} furthest points from vertex 1 to vertex 1, second \code{k} entries are distances from the \code{k} furthest
points from vertex 2 to vertex 2, and the last \code{k} entries are the distances from the \code{k} furthest points
from vertex 3 to vertex 3.}
}
\description{
An object of class "Extrema".
Returns the \code{k} furthest data points among the data set, \code{Dt}, in each CC-vertex region from the vertex in the
triangle, \eqn{tri=T(A,B,C)}, vertices are stacked row-wise. Vertex region labels/numbers correspond to the
row number of the vertex in \code{tri}.

\code{ch.all.intri} is for checking whether all data points are inside \code{tri} (default is \code{FALSE}).
If some of the data points are not inside \code{tri} and \code{ch.all.intri}=\code{TRUE}, then the function yields
an error message.
If some of the data points are not inside \code{tri} and \code{ch.all.intri}=\code{FALSE}, then the function yields
the closest points to edges among the data points inside \code{tri} (yields \code{NA} if there are no data points
inside \code{tri}).

In the extrema, \eqn{ext}, in the output, the first \code{k} entries are the \code{k} furthest points from vertex 1,
second \code{k} entries are \code{k} furthest points are from vertex 2, and
last \code{k} entries are the \code{k} furthest points from vertex 3.
If data size does not allow, NA's are inserted for some or all of the \code{k} furthest points for each vertex.
}
\examples{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10  #try also n<-20
k<-3

set.seed(1)
dat<-runif.tri(n,Tr)$g

Ext<-Kfr2vVRCC(dat,Tr,k)
Ext
summary(Ext)
plot(Ext)

Kfr2vVRCC(dat[1:k,],Tr,k)
Kfr2vVRCC(dat[1,],Tr,k)

dat2<-rbind(dat,c(.2,.4))
Kfr2vVRCC(dat2,Tr,k)  #try also Kfr2vVRCC(dat2,Tr,k,ch.all.intri = TRUE)

kf2v<-Kfr2vVRCC(dat,Tr,k)

CC<-circ.cent.tri(Tr)  #the circumcenter
D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
Ds<-rbind(D1,D2,D3)

Xlim<-range(Tr[,1],dat[,1])
Ylim<-range(Tr[,2],dat[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(A,pch=".",asp=1,xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty=2)
points(dat)
points(kf2v$ext,pch=4,col=2)

txt<-rbind(Tr,CC,Ds)
xc<-txt[,1]+c(-.06,.08,.05,.12,-.1,-.1,-.09)
yc<-txt[,2]+c(.02,-.02,.04,.0,.02,.06,-.04)
txt.str<-c("A","B","C","CC","D1","D2","D3")
text(xc,yc,txt.str)

Kfr2vVRCC(dat,Tr,k)

dat.fr<-data.frame(a=dat)
Kfr2vVRCC(dat.fr,Tr,k)

dat.fr<-data.frame(a=Tr)
fr2vVRCC(dat,dat.fr)

}
\seealso{
\code{\link{Kfr2vTbVRCC}}, \code{\link{fr2vTbVRCC}}, \code{\link{fr2vVRCC}}, and \code{\link{fr2eTeER}}
}
