% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PropEdge3D.R
\name{Idom.num1PEtetra}
\alias{Idom.num1PEtetra}
\title{The indicator for a 3D point being a dominating point for Proportional Edge Proximity Catch Digraphs (PE-PCDs)
- one tetrahedron case}
\usage{
Idom.num1PEtetra(p, Xp, th, r, M = "CM", rv = NULL, ch.data.pnt = FALSE)
}
\arguments{
\item{p}{A 3D point that is to be tested for being a dominating point or not of the PE-PCD.}

\item{Xp}{A set of 3D points which constitutes the vertices of the PE-PCD.}

\item{th}{A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.}

\item{r}{A positive real number which serves as the expansion parameter in PE proximity region;
must be \eqn{\ge 1}.}

\item{M}{The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.}

\item{rv}{Index of the vertex whose region contains point \code{p}, \code{rv} takes the vertex labels as \code{1,2,3,4} as
in the row order of the vertices in standard tetrahedron, default is \code{NULL}.}

\item{ch.data.pnt}{A logical argument for checking whether point \code{p} is a data point
in \code{Xp} or not (default is \code{FALSE}).}
}
\value{
\eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 2D data set \code{Xp},
that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
}
\description{
Returns \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 2D data set \code{Xp} in the tetrahedron \code{th}, that is,
returns 1 if \code{p} is a dominating point of PE-PCD, returns 0 otherwise.

Point, \code{p}, is in the vertex region of vertex \code{rv} (default is \code{NULL}); vertices are labeled as \code{1,2,3,4}
in the order they are stacked row-wise in \code{th}.

PE proximity region is constructed with respect to the tetrahedron \code{th} with expansion parameter \eqn{r \ge 1}
and vertex regions are based on center of mass (\code{M="CM"}) or circumcenter (\code{M="CC"}) only.
and vertex regions are based on center of mass \eqn{CM} (equivalent to circumcenter in this case).

\code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
so by default this function checks whether the point \code{p} would be a dominating point
if it actually were in the data set.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
tetra<-rbind(A,B,C,D)
n<-5 #try also n<-20

Xp<-runif.tetra(n,tetra)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))

M<-"CM"; cent<-apply(tetra,2,mean)  #center of mass
#try also M<-"CC"; cent<-circumcenter.tetra(tetra)  #circumcenter

r<-2

P<-c(.4,.1,.2)
Idom.num1PEtetra(Xp[1,],Xp,tetra,r,M)
Idom.num1PEtetra(P,Xp,tetra,r,M)

#or try
RV<-rel.vert.tetraCC(Xp[1,],tetra)$rv
Idom.num1PEtetra(Xp[1,],Xp,tetra,r,M,rv=RV)

Idom.num1PEtetra(c(-1,-1,-1),Xp,tetra,r,M)
Idom.num1PEtetra(c(-1,-1,-1),c(-1,-1,-1),tetra,r,M)

gam.vec<-vector()
for (i in 1:n)
{gam.vec<-c(gam.vec,Idom.num1PEtetra(Xp[i,],Xp,tetra,r,M))}

ind.gam1<-which(gam.vec==1)
ind.gam1
g1.pts<-Xp[ind.gam1,]

Xlim<-range(tetra[,1],Xp[,1],cent[1])
Ylim<-range(tetra[,2],Xp[,2],cent[2])
Zlim<-range(tetra[,3],Xp[,3],cent[3])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]
zd<-Zlim[2]-Zlim[1]

plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi =0,theta=40, bty = "g",
xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
         pch = 20, cex = 1, ticktype = "detailed")
#add the vertices of the tetrahedron
plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
if (length(g1.pts)!=0)
{plot3D::points3D(g1.pts[,1],g1.pts[,2],g1.pts[,3], pch=4,col="red", add=TRUE)}

plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)

D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
L<-rbind(D1,D2,D3,D4,D5,D6); R<-rbind(cent,cent,cent,cent,cent,cent)
plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)

P<-c(.4,.1,.2)
Idom.num1PEtetra(P,Xp,tetra,r,M)

Idom.num1PEtetra(c(-1,-1,-1),Xp,tetra,r,M,ch.data.pnt = FALSE)
#gives an error message if ch.data.pnt = TRUE since p is not a data point
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{Idom.num1PEstd.tetra}}, \code{\link{Idom.num1PEtri}} and \code{\link{Idom.num1PEbasic.tri}}
}
\author{
Elvan Ceyhan
}
