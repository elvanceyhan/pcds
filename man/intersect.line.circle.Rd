% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceFunctions.R
\name{intersect.line.circle}
\alias{intersect.line.circle}
\title{The points of intersection of a line and a circle}
\usage{
intersect.line.circle(p1, p2, cent, rad)
}
\arguments{
\item{p1, p2}{2D points that determine the straight line (i.e., through which the straight line passes).}

\item{cent}{A 2D point representing the center of the circle.}

\item{rad}{A positive real number representing the radius of the circle.}
}
\value{
point(s) of intersection between the circle and the line (if they do not intersect, the function
yields \code{NULL} as the output)
}
\description{
Returns the intersection point(s) of a line and a circle. The line is determined by the two points
\code{p1} and \code{p2} and the circle is centered at point \code{cent} and has radius \code{rad}.
If the circle does not intersect the line, the function yields \code{NULL};
if the circle intersects at only one point, it yields only that point;
otherwise it yields both intersection points as output.
When there are two intersection points, they are listed in the order of the \eqn{x}-coordinates of \code{p1} and \code{p2};
and if the \eqn{x}-coordinates of \code{p1} and \code{p2} are equal, intersection points are listed in the order of
\eqn{y}-coordinates of \code{p1} and \code{p2}.
}
\examples{
\dontrun{
P1<-c(.3,.2)*100
P2<-c(.6,.3)*100
cent<-c(1.1,1.1)*100
rad<-2*100

intersect.line.circle(P1,P2,cent,rad)
intersect.line.circle(P2,P1,cent,rad)
intersect.line.circle(P1,P1+c(0,1),cent,rad)
intersect.line.circle(P1+c(0,1),P1,cent,rad)

dist.point2line(cent,P1,P2)
rad2<-dist.point2line(cent,P1,P2)$d
intersect.line.circle(P1,P2,cent,rad2)
intersect.line.circle(P1,P2,cent,rad=.8)
intersect.line.circle(P1,P2,cent,rad=.78)

#plot of the line and the circle
A<-c(.3,.2); B<-c(.6,.3); cent<-c(1,1); rad<-2 #check dist.point2line(cent,A,B)$dis, .3

IPs<-intersect.line.circle(A,B,cent,rad)

xr<-range(A[1],B[1],cent[1])
xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
x<-seq(xr[1]-rad-xf,xr[2]+rad+xf,l=20)  #try also l=100
lnAB<-Line(A,B,x)
y<-lnAB$y

Xlim<-range(x,cent[1])
Ylim<-range(y,A[2],B[2],cent[2]-rad,cent[2]+rad)
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(rbind(A,B,cent),pch=1,asp=1,xlab="x",ylab="y",
xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
lines(x,y,lty=1)
interp::circles(cent[1],cent[2],rad)
IP.txt<-c()
if (!is.null(IPs))
{
  for (i in 1:(length(IPs)/2))
    IP.txt<-c(IP.txt,paste("I",i, sep = ""))
}
txt<-rbind(A,B,cent,IPs)
text(txt+cbind(rep(xd*.03,nrow(txt)),rep(-yd*.03,nrow(txt))),c("A","B","M",IP.txt))
}

}
\seealso{
\code{\link{intersect2lines}}
}
\author{
Elvan Ceyhan
}
