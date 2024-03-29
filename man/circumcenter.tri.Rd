% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxDelaunay.R
\name{circumcenter.tri}
\alias{circumcenter.tri}
\title{Circumcenter of a general triangle}
\usage{
circumcenter.tri(tri)
}
\arguments{
\item{tri}{A \eqn{3 \times 2} matrix with each row
representing a vertex of the triangle.}
}
\value{
circumcenter of the triangle \code{tri}
}
\description{
Returns the circumcenter a given triangle, \code{tri},
with vertices stacked row-wise.
See (\insertCite{weisstein-tri-centers,ceyhan:comp-geo-2010;textual}{pcds})
for triangle centers.
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);  #the vertices of the triangle Tr

CC<-circumcenter.tri(Tr)  #the circumcenter
CC

D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
Ds<-rbind(D1,D2,D3)

Xlim<-range(Tr[,1],CC[1])
Ylim<-range(Tr[,2],CC[2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(A,asp=1,pch=".",xlab="",ylab="",main="Circumcenter of a triangle",
axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
points(rbind(CC))
L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty = 2)

txt<-rbind(Tr,CC,Ds)
xc<-txt[,1]+c(-.08,.08,.08,.12,-.09,-.1,-.09)
yc<-txt[,2]+c(.02,-.02,.03,-.06,.02,.06,-.04)
txt.str<-c("A","B","C","CC","D1","D2","D3")
text(xc,yc,txt.str)

A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
Te<-rbind(A,B,C);  #the vertices of the equilateral triangle Te
circumcenter.tri(Te)  #the circumcenter

A<-c(0,0); B<-c(0,1); C<-c(2,0);
Tr<-rbind(A,B,C);  #the vertices of the triangle T
circumcenter.tri(Tr)  #the circumcenter
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{circumcenter.basic.tri}}
}
\author{
Elvan Ceyhan
}
