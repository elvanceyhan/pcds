% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxDelaunay.R
\name{rel.verts.triM}
\alias{rel.verts.triM}
\title{The alternative function for the indices of the M-vertex regions
in a triangle that contains the points
in a give data set}
\usage{
rel.verts.triM(Xp, tri, M)
}
\arguments{
\item{Xp}{A set of 2D points representing the set of data points
for which indices of the vertex regions
containing them are to be determined.}

\item{tri}{A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.}

\item{M}{A 2D point in Cartesian coordinates
or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri}.}
}
\value{
A \code{list} with two elements
\item{rv}{Indices of the vertices
whose regions contains points in \code{Xp}.}
\item{tri}{The vertices of the triangle,
where row number corresponds to the vertex index in \code{rv}.}
}
\description{
An alternative function to the function \code{\link{rel.verts.tri}}
when the center M is not the circumcenter falling outside the triangle.
This function only works for a center \eqn{M} in the interior of the triangle,
with the projections of \eqn{M} to the edges
along the lines joining \eqn{M} to the vertices.
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
M<-c(1.6,1.0)

P<-c(.4,.2)
rel.verts.triM(P,Tr,M)

n<-20  #try also n<-40
set.seed(1)
Xp<-runif.tri(n,Tr)$g

M<-c(1.6,1.0)  #try also M<-c(1.3,1.3)

(rv<-rel.verts.tri(Xp,Tr,M))
rel.verts.triM(rbind(Xp,c(2,2)),Tr,M)

Ds<-prj.cent2edges(Tr,M)

Xlim<-range(Tr[,1])
Ylim<-range(Tr[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
points(Xp,pch=".",col=1)
L<-rbind(M,M,M); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty = 2)

xc<-Tr[,1]+c(-.03,.05,.05)
yc<-Tr[,2]+c(-.06,.02,.05)
txt.str<-c("rv=1","rv=2","rv=3")
text(xc,yc,txt.str)

txt<-rbind(M,Ds)
xc<-txt[,1]+c(.02,.04,-.03,0)
yc<-txt[,2]+c(.07,.03,.05,-.07)
txt.str<-c("M","D1","D2","D3")
text(xc,yc,txt.str)
text(Xp,labels=factor(rv$rv))
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{rel.verts.tri}}
}
\author{
Elvan Ceyhan
}
