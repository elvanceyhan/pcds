% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxGeometry.R
\name{paraplane}
\alias{paraplane}
\title{The plane at a point and parallel to the plane spanned
by three distinct 3D points \code{a}, \code{b}, and \code{c}}
\usage{
paraplane(p, a, b, c, x, y)
}
\arguments{
\item{p}{A 3D point which the plane parallel to the plane spanned by
three distinct 3D points \code{a}, \code{b}, and \code{c} crosses.}

\item{a, b, c}{3D points that determine the plane to which the plane
crossing point \code{p} is parallel to.}

\item{x, y}{Scalars or \code{vectors} of scalars
representing the \eqn{x}- and \eqn{y}-coordinates of the plane parallel to
the plane spanned by points \code{a}, \code{b},
and \code{c} and passing through point \code{p}.}
}
\value{
A \code{list} with the elements
\item{desc}{Description of the plane passing through point \code{p}
and parallel to plane spanned by points
\code{a}, \code{b} and \code{c}}
\item{points}{The input points \code{a}, \code{b}, \code{c}, and \code{p}.
Plane is parallel to the plane spanned by \code{a}, \code{b}, and \code{c}
and passes through point \code{p} (stacked row-wise,
i.e., row 1 is point \code{a}, row 2 is point \code{b},
 row 3 is point \code{c}, and row 4 is point \code{p}).}
\item{x,y}{The input vectors which constitutes the \eqn{x}-
and \eqn{y}-coordinates of the point(s) of interest on the
plane. \code{x} and \code{y} can be scalars or vectors of scalars.}
\item{z}{The output \code{vector}
which constitutes the \eqn{z}-coordinates of the point(s) of interest on the plane.
If \code{x} and \code{y} are scalars, \code{z} will be a scalar and
if \code{x} and \code{y} are vectors of scalars,
then \code{z} needs to be a \code{matrix} of scalars,
containing the \eqn{z}-coordinate for each pair of \code{x} and \code{y} values.}
\item{coeff}{Coefficients of the plane (in the \eqn{z = A x+B y+C} form).}
\item{equation}{Equation of the plane in long form}
\item{equation2}{Equation of the plane in short form,
to be inserted on the plot}
}
\description{
An object of class \code{"Planes"}.
Returns the equation and \eqn{z}-coordinates of the plane
passing through point \code{p} and parallel to the plane spanned
by three distinct 3D points \code{a}, \code{b},
and \code{c} with \eqn{x}- and \eqn{y}-coordinates are provided
in vectors \code{x} and \code{y},
respectively.
}
\examples{
\dontrun{
Q<-c(1,10,3); R<-c(1,1,3); S<-c(3,9,12); P<-c(1,1,0)

pts<-rbind(Q,R,S,P)
paraplane(P,Q,R,S,.1,.2)

xr<-range(pts[,1]); yr<-range(pts[,2])
xf<-(xr[2]-xr[1])*.25
#how far to go at the lower and upper ends in the x-coordinate
yf<-(yr[2]-yr[1])*.25
#how far to go at the lower and upper ends in the y-coordinate
x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
y<-seq(yr[1]-yf,yr[2]+yf,l=5)  #try also l=10, 20, or 100

plP2QRS<-paraplane(P,Q,R,S,x,y)
plP2QRS
summary(plP2QRS)
plot(plP2QRS,theta = 225, phi = 30, expand = 0.7, facets = FALSE, scale = TRUE)

paraplane(P,Q,R,Q+R,.1,.2)

z.grid<-plP2QRS$z

plQRS<-Plane(Q,R,S,x,y)
plQRS
pl.grid<-plQRS$z

zr<-max(z.grid)-min(z.grid)
Pts<-rbind(Q,R,S,P)+rbind(c(0,0,zr*.1),c(0,0,zr*.1),
c(0,0,zr*.1),c(0,0,zr*.1))
Mn.pts<-apply(Pts[1:3,],2,mean)

plot3D::persp3D(z = pl.grid, x = x, y = y, theta =225, phi = 30,
ticktype = "detailed",
main="Plane Crossing Points Q, R, S\n and Plane Passing P Parallel to it")
#plane spanned by points Q, R, S
plot3D::persp3D(z = z.grid, x = x, y = y,add=TRUE)
#plane parallel to the original plane and passing thru point \code{P}

plot3D::persp3D(z = z.grid, x = x, y = y, theta =225, phi = 30,
ticktype = "detailed",
main="Plane Crossing Point P \n and Parallel to the Plane Crossing Q, R, S")
#plane spanned by points Q, R, S
#add the defining points
plot3D::points3D(Pts[,1],Pts[,2],Pts[,3], add=TRUE)
plot3D::text3D(Pts[,1],Pts[,2],Pts[,3], c("Q","R","S","P"),add=TRUE)
plot3D::text3D(Mn.pts[1],Mn.pts[2],Mn.pts[3],plP2QRS$equation,add=TRUE)
plot3D::polygon3D(Pts[1:3,1],Pts[1:3,2],Pts[1:3,3], add=TRUE)
}

}
\seealso{
\code{\link{Plane}}
}
\author{
Elvan Ceyhan
}
