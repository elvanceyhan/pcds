% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceFunctions.R
\name{Idom.num1AStri}
\alias{Idom.num1AStri}
\title{The indicator for a point being a dominating point for Arc Slice Proximity Catch Digraphs
(AS-PCDs) - one triangle case}
\usage{
Idom.num1AStri(p, Xp, tri, M = "CC", rv = NULL, ch.data.pnt = FALSE)
}
\arguments{
\item{p}{A 2D point that is to be tested for being a dominating point or not of the AS-PCD.}

\item{Xp}{A set of 2D points which constitutes the vertices of the AS-PCD.}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{M}{The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
a 3D point in  barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
default is \code{M="CC"} i.e., the circumcenter of \code{tri}.}

\item{rv}{Index of the vertex whose region contains point \code{p}, \code{rv} takes the vertex labels as \eqn{1,2,3} as
in the row order of the vertices in \code{tri}.}

\item{ch.data.pnt}{A logical argument for checking whether point \code{p} is a data point in \code{Xp} or not
(default is \code{FALSE}).}
}
\value{
I(\code{p} is a dominating point of the AS-PCD whose vertices are the 2D data set \code{Xp}),
that is, returns 1 if \code{p} is a dominating point of the AS-PCD, returns 0 otherwise
}
\description{
Returns I(\code{p} is a dominating point of the AS-PCD whose vertices are the 2D data set \code{Xp}),
that is, returns 1 if \code{p} is a dominating point of AS-PCD, returns 0 otherwise.
Point, \code{p}, is in the region of vertex \code{rv} (default is \code{NULL}); vertices are labeled as \eqn{1,2,3}
in the order they are stacked row-wise in \code{tri}.

AS proximity regions are defined with respect to the
triangle \code{tri} and vertex regions are based on the center \code{M="CC"} for circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
default is \code{M="CC"} the circumcenter of \code{tri}.

\code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
so by default this function checks whether the point \code{p} would be a dominating point
if it actually were in the data set.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10

set.seed(1)
Xp<-runif.tri(n,Tr)$g

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)

Idom.num1AStri(Xp[1,],Xp,Tr,M)
Idom.num1AStri(Xp[1,],Xp[1,],Tr,M)
Idom.num1AStri(c(1.5,1.5),c(1.6,1),Tr,M)
Idom.num1AStri(c(1.6,1),c(1.5,1.5),Tr,M)

gam.vec<-vector()
for (i in 1:n)
{gam.vec<-c(gam.vec,Idom.num1AStri(Xp[i,],Xp,Tr,M))}

ind.gam1<-which(gam.vec==1)
ind.gam1

#or try
Rv<-rel.vert.triCC(Xp[1,],Tr)$rv
Idom.num1AStri(Xp[1,],Xp,Tr,M,Rv)

Idom.num1AStri(c(.2,.4),Xp,Tr,M)
Idom.num1AStri(c(.2,.4),c(.2,.4),Tr,M)

Xp2<-rbind(Xp,c(.2,.4))
Idom.num1AStri(Xp[1,],Xp2,Tr,M)

if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#need to run this when M is given in barycentric coordinates

CC<-circumcenter.tri(Tr)  #the circumcenter

if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
{cent<-CC
D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
Ds<-rbind(D1,D2,D3)
cent.name<-"CC"
} else
{cent<-M
cent.name<-"M"
Ds<-prj.cent2edges(Tr,M)
}

Xlim<-range(Tr[,1],Xp[,1])
Ylim<-range(Tr[,2],Xp[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(A,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
points(Xp)
L<-rbind(cent,cent,cent); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty=2)
points(rbind(Xp[ind.gam1,]),pch=4,col=2)

txt<-rbind(Tr,cent,Ds)
xc<-txt[,1]
yc<-txt[,2]
txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
text(xc,yc,txt.str)

Idom.num1AStri(c(1.5,1.1),Xp,Tr,M)

Idom.num1AStri(c(1.5,1.1),Xp,Tr,M)

Idom.num1AStri(c(1.5,1.1),Xp,Tr,M,ch.data.pnt=FALSE)
#gives an error message if ch.data.pnt=TRUE since point p is not a data point in Xp
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{Idom.num1ASbasic.tri}}
}
\author{
Elvan Ceyhan
}
