% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxRFuncs4PCDClasses.r
\name{print.Planes}
\alias{print.Planes}
\title{Print a \code{Planes} object}
\usage{
\method{print}{Planes}(x, ...)
}
\arguments{
\item{x}{A \code{Planes} object.}

\item{\dots}{Additional arguments for the S3 method '\code{print}'.}
}
\value{
The \code{call} of the object of class '\code{Planes}'
and the \code{coefficients} of the plane (in the form: z = A*x + B*y + C).
}
\description{
Prints the \code{call} of the object of class '\code{Planes}'
and also the \code{coefficients} of the plane (in the form: z = A*x + B*y + C).
}
\examples{
A<-c(1,10,3); B<-c(1,1,3); C<-c(3,9,12)
pts<-rbind(A,B,C)

xr<-range(pts[,1]); yr<-range(pts[,2])
xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
yf<-(yr[2]-yr[1])*.1 #how far to go at the lower and upper ends in the y-coordinate
x<-seq(xr[1]-xf,xr[2]+xf,l=100)
y<-seq(yr[1]-yf,yr[2]+yf,l=100)

plABC<-Plane(A,B,C,x,y)
plABC
print(plABC)

}
\seealso{
\code{\link{summary.Planes}}, \code{\link{print.summary.Planes}}, and \code{\link{plot.Planes}}
}
