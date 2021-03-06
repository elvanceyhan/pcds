% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{cl2Mc.int}
\alias{cl2Mc.int}
\title{The closest points to center in each vertex region in an interval}
\usage{
cl2Mc.int(Dt, int, c)
}
\arguments{
\item{Dt}{A set or vector of 1D points from which closest points to \eqn{M_c} are found
in the interval \code{int}.}

\item{int}{A vector of two real numbers representing an interval.}

\item{c}{A positive real number in \eqn{(0,1)} parameterizing the center inside \eqn{int=(a,b)}.
For the interval, \eqn{int=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.}
}
\value{
A list with the elements
\item{txt1}{Vertex Labels are a=1 and b=2 for the interval (a,b).}
\item{txt2}{A short description of the distances as "Distances from ..."}
\item{type}{Type of the extrema points}
\item{desc}{A short description of the extrema points}
\item{mtitle}{The "main" title for the plot of the exrema}
\item{ext}{The extrema points, here, closest points to \eqn{M_c} in each vertex region}
\item{X}{The input data vector, \code{Dt}.}
\item{num.points}{The number of data points, i.e., size of \code{Dt}}
\item{supp}{Support of the data points, here, it is int.}
\item{cent}{The (parameterized) center point used for construction of vertex regions.}
\item{ncent}{Name of the (parameterized) center, \code{cent}, it is "Mc" for this function.}
\item{regions}{Vertex regions inside the interval, \code{int}, provided as a list.}
\item{region.names}{Names of the vertex regions as "vr=1","vr=2"}
\item{region.centers}{Centers of mass of the vertex regions inside \code{int}.}
\item{dist2ref}{Distances from closest points in each vertex region to \eqn{M_c}.}
}
\description{
An object of class "Extrema".
Returns the closest data points among the data set, \code{Dt}, in each \eqn{M_c}-vertex region
i.e., finds the closest points from right and left to \eqn{M_c} among points of the 1D data set \code{Dt} which reside in
in the interval \eqn{int=(a,b)}.

\eqn{M_c} is based on the centrality parameter \code{c} in \eqn{(0,1)},
so that \eqn{100c} \% of the length of interval is to the left of \eqn{M_c} and \eqn{100(1-c)} \% of the length of the interval
is to the right of \eqn{M_c}. That is, for the interval \eqn{(a,b)}, \eqn{M_c=a+c(b-a)}. If there are no points from Dt to
the left of \eqn{M_c} in the interval, then it yields NA, and likewise for the right of \eqn{M_c} in the interval.

See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
}
\examples{
c<-.4
a<-0; b<-10; int<-c(a,b)

Mc<-centMc(int,c)

nx<-10
xr<-range(a,b,Mc)
xf<-(xr[2]-xr[1])*.5

dat<-runif(nx,a,b)

Ext<-cl2Mc.int(dat,int,c)
Ext
summary(Ext)
plot(Ext)

cl2Mc.int(dat[1],int,c)
cMc<-cl2Mc.int(dat,int,c)

Xlim<-range(a,b,dat)
xd<-Xlim[2]-Xlim[1]

plot(cbind(a,0),xlab="",pch=".",xlim=Xlim+xd*c(-.05,.05))
abline(h=0)
abline(v=c(a,b,Mc),col=c(1,1,2),lty=2)
points(cbind(dat,0))
points(cbind(c(cMc$ext),0),pch=4,col=2)
text(cbind(c(a,b,Mc),-0.1),c("a","b","Mc"))

n<-10  #try also n<-20
dat<-runif(n,a-5,b+5)
cl2Mc.int(dat,c(a,b),c)

dat<-runif(n,a+b,b+10)
cl2Mc.int(dat,int,c)

c<-.4
a<-0; b<-10; int<-c(a,b)
n<-10
dat<-runif(n,a,b)
cl2Mc.int(dat,int,c)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{cl2CC.TbVR}} and \code{\link{cl2CC.VR}}
}
