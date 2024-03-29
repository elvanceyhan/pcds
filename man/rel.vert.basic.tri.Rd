% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxDelaunay.R
\name{rel.vert.basic.tri}
\alias{rel.vert.basic.tri}
\title{The index of the vertex region in a standard basic triangle form
that contains a given point}
\usage{
rel.vert.basic.tri(p, c1, c2, M)
}
\arguments{
\item{p}{A 2D point for which \code{M}-vertex region it resides in is
to be determined in the
standard basic triangle form \eqn{T_b}.}

\item{c1, c2}{Positive real numbers
which constitute the vertex of the standard basic triangle form
adjacent to the shorter edges;
\eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0}
and \eqn{(1-c_1)^2+c_2^2 \le 1}.}

\item{M}{A 2D point in Cartesian coordinates
or a 3D point in barycentric coordinates
which serves as a center
in the interior of the standard basic triangle form.}
}
\value{
A \code{list} with two elements
\item{rv}{Index of the vertex whose region contains point, \code{p};
index of the vertex is the same as the row
number in the standard basic triangle form, \eqn{T_b}}
\item{tri}{The vertices of the standard basic triangle form, \eqn{T_b},
where row number corresponds to the vertex index \code{rv}
with \code{rv=1} is row \eqn{1=(0,0)}, \code{rv=2} is row \eqn{2=(1,0)},
and \eqn{rv=3} is row \eqn{3=(c_1,c_2)}.}
}
\description{
Returns the index of the related vertex
in the standard basic triangle form
whose region contains point \code{p}.
The standard basic triangle form is \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}..

Vertex regions are based on the general center \eqn{M=(m_1,m_2)}
in Cartesian coordinates or
\eqn{M=(\alpha,\beta,\gamma)}
in barycentric coordinates
in the interior of the standard basic triangle form \eqn{T_b}.
Vertices of the standard basic triangle form \eqn{T_b} are labeled
according to the row number the
vertex is recorded, i.e., as 1=(0,0), 2=(1,0),and \eqn{3=(c_1,c_2)}.

If the point, \code{p}, is not inside \eqn{T_b},
then the function yields \code{NA} as output.
The corresponding vertex region is the polygon
with the vertex, \code{M}, and projections from \code{M}
to the edges on the lines joining vertices and \code{M}.
That is, \code{rv=1} has vertices \eqn{(0,0),D_3,M,D_2};
\code{rv=2} has vertices \eqn{(1,0),D_1,M,D_3}; and
\eqn{rv=3} has vertices \eqn{(c_1,c_2),D_2,M,D_1}
(see the illustration in the examples).

Any given triangle can be mapped to the standard basic triangle form
by a combination of rigid body motions
(i.e., translation, rotation and reflection) and scaling,
preserving uniformity of the points in the original triangle.
Hence, standard basic triangle form is useful for simulation
studies under the uniformity hypothesis.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
}
\examples{
\dontrun{
c1<-.4; c2<-.6
A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
Tb<-rbind(A,B,C);
M<-c(.6,.2)

P<-c(.4,.2)
rel.vert.basic.tri(P,c1,c2,M)

n<-20  #try also n<-40
set.seed(1)
Xp<-runif.basic.tri(n,c1,c2)$g

M<-as.numeric(runif.basic.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)

Rv<-vector()
for (i in 1:n)
{ Rv<-c(Rv,rel.vert.basic.tri(Xp[i,],c1,c2,M)$rv)}
Rv

Ds<-prj.cent2edges.basic.tri(c1,c2,M)

Xlim<-range(Tb[,1],Xp[,1])
Ylim<-range(Tb[,2],Xp[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

if (dimension(M)==3) {M<-bary2cart(M,Tb)}
#need to run this when M is given in barycentric coordinates

plot(Tb,pch=".",xlab="",ylab="",axes=TRUE,
xlim=Xlim+xd*c(-.1,.1),ylim=Ylim+yd*c(-.05,.05))
polygon(Tb)
points(Xp,pch=".",col=1)
L<-rbind(M,M,M); R<-Ds
segments(L[,1], L[,2], R[,1], R[,2], lty = 2)

xc<-Tb[,1]+c(-.04,.05,.04)
yc<-Tb[,2]+c(.02,.02,.03)
txt.str<-c("rv=1","rv=2","rv=3")
text(xc,yc,txt.str)

txt<-rbind(M,Ds)
xc<-txt[,1]+c(-.02,.04,-.03,0)
yc<-txt[,2]+c(-.02,.02,.02,-.03)
txt.str<-c("M","D1","D2","D3")
text(xc,yc,txt.str)
text(Xp,labels=factor(Rv))
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{rel.vert.basic.triCM}}, \code{\link{rel.vert.tri}},
\code{\link{rel.vert.triCC}}, \code{\link{rel.vert.basic.triCC}},
\code{\link{rel.vert.triCM}}, and \code{\link{rel.vert.std.triCM}}
}
\author{
Elvan Ceyhan
}
