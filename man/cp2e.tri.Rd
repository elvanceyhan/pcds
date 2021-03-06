% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{cp2e.tri}
\alias{cp2e.tri}
\title{Projections of a point inside a triangle to its edges}
\usage{
cp2e.tri(tri, M)
}
\arguments{
\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri}.}
}
\value{
Three projection points (stacked row-wise) from a general center \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of a triangle to the edges on
the extension of the lines joining \code{M} to the vertices;
row \eqn{i} is the projection  point into edge \eqn{i}, for \eqn{i=1,2,3}.
}
\description{
Returns the projections from a general center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
barycentric coordinates in the interior of a triangle to the edges on the extension of the lines joining
\code{M} to the vertices (see the examples for an illustration).

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)

cp2e.tri(Tr,M)  #try also cp2e.tri(Tr,M=c(1,1))

Ds<-cp2e.tri(Tr,M)

Xlim<-range(Tr[,1])
Ylim<-range(Tr[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#need to run this when M is given in barycentric coordinates

plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,
xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
L<-rbind(M,M,M); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty=2)

xc<-Tr[,1]
yc<-Tr[,2]
txt.str<-c("rv=1","rv=2","rv=3")
text(xc,yc,txt.str)

txt<-rbind(M,Ds)
xc<-txt[,1]+c(-.02,.04,-.04,-.02)
yc<-txt[,2]+c(-.02,.04,.04,-.06)
txt.str<-c("M","D1","D2","D3")
text(xc,yc,txt.str)

cp2e.tri(Tr,M)

dat.fr<-data.frame(a=Tr)
cp2e.tri(dat.fr,M)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{cp2e.bastri}} and \code{\link{cp2edges.nd}}
}
