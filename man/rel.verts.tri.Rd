% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxDelaunay.R
\name{rel.verts.tri}
\alias{rel.verts.tri}
\title{The indices of the vertex regions in a triangle
that contains the points in a give data set}
\usage{
rel.verts.tri(Xp, tri, M)
}
\arguments{
\item{Xp}{A set of 2D points representing the set of data points
for which indices of the vertex regions containing them are to be determined.}

\item{tri}{A \eqn{3 \times 2} matrix
with each row representing a vertex of the triangle.}

\item{M}{A 2D point in Cartesian coordinates
or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri}
or the circumcenter of \code{tri}.}
}
\value{
A \code{list} with two elements
\item{rv}{Indices of the vertices
whose regions contains points in \code{Xp}.}
\item{tri}{The vertices of the triangle,
where row number corresponds to the vertex index in \code{rv}.}
}
\description{
Returns the indices of the vertices
whose regions contain the points in data set \code{Xp} in
a triangle \code{tri}\eqn{=T(A,B,C)}.

Vertex regions are based on center \eqn{M=(m_1,m_2)}
in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
in barycentric coordinates in the interior of the triangle
to the edges on the extension of the lines joining \code{M} to the vertices
or based on the circumcenter of \code{tri}.
Vertices of triangle \code{tri} are labeled as \eqn{1,2,3}
according to the row number the vertex is recorded.

If a point in \code{Xp} is not inside \code{tri},
then the function yields \code{NA} as output for that entry.
The corresponding vertex region is the polygon with the vertex, \code{M}, and
projection points from \code{M} to the edges crossing the vertex
(as the output of \code{prj.cent2edges(Tr,M)})
or \eqn{CC}-vertex region
(see the examples for an illustration).

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:dom-num-NPE-Spat2011,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
M<-c(1.6,1.0)

P<-c(.4,.2)
rel.verts.tri(P,Tr,M)

n<-20  #try also n<-40
set.seed(1)
Xp<-runif.tri(n,Tr)$g

M<-as.numeric(runif.tri(1,Tr)$g)  #try also  #M<-c(1.6,1.0)

rel.verts.tri(Xp,Tr,M)
rel.verts.tri(rbind(Xp,c(2,2)),Tr,M)

rv<-rel.verts.tri(Xp,Tr,M)
rv

ifelse(identical(M,circumcenter.tri(Tr)),
Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2),Ds<-prj.cent2edges(Tr,M))

Xlim<-range(Tr[,1],M[1],Xp[,1])
Ylim<-range(Tr[,2],M[2],Xp[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#need to run this when M is given in barycentric coordinates

plot(Tr,pch=".",xlab="",ylab="",
main="Scatterplot of data points \n and M-vertex regions in a triangle",
axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
points(Xp,pch=".",col=1)
L<-rbind(M,M,M); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty = 2)

xc<-Tr[,1]
yc<-Tr[,2]
txt.str<-c("rv=1","rv=2","rv=3")
text(xc,yc,txt.str)

txt<-rbind(M,Ds)
xc<-txt[,1]+c(.02,.04,-.03,0)
yc<-txt[,2]+c(.07,.04,.05,-.07)
txt.str<-c("M","D1","D2","D3")
text(xc,yc,txt.str)
text(Xp,labels=factor(rv$rv))
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{rel.verts.triCM}}, \code{\link{rel.verts.triCC}},
and \code{\link{rel.verts.tri.nondegPE}}
}
\author{
Elvan Ceyhan
}
