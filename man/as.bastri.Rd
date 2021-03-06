% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{as.bastri}
\alias{as.bastri}
\title{The labels of the vertices of a triangle in the (unscaled) basic triangle form}
\usage{
as.bastri(tri)
}
\arguments{
\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}
}
\value{
A list with two elements
\item{tri}{The vertices of the triangle stacked row-wise and labeled row-wise as A, B, C.}
\item{desc}{Description of the edges based on the vertices, i.e.,
 "Edges (in decreasing length are) AB, BC, and AC".}
\item{orig.order}{Row order of the input triangle, \code{tri},
when converted to the scaled version of the basic triangle}
}
\description{
Labels the vertices of triangle, \code{tri}, as ABC so that AB is the longest edge, BC is the second longest
and AC is the shortest edge (the order is as in the basic triangle). The new triangle \eqn{T(A,B,C)} is unscaled,
i.e., the longest edge AB may not be of unit length.

The basic triangle is \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} where \code{c1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}

Any given triangle can be mapped to the basic triangle by a combination of rigid body motions
(i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
original triangle. Hence basic triangle is useful for simulation
studies under the uniformness hypothesis.
}
\examples{
\dontrun{
c1<-.4; c2<-.6
A<-c(0,0); B<-c(1,0); C<-c(c1,c2);

as.bastri(rbind(A,B,C))

as.bastri(rbind(B,C,A))

as.bastri(rbind(B,A,C))
as.bastri(rbind(A,C,B))

A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
as.bastri(rbind(A,B,C))
as.bastri(rbind(A,C,B))
as.bastri(rbind(B,A,C))

A<-runif(2); B<-runif(2); C<-runif(2)
as.bastri(rbind(A,B,C))

dat.fr<-data.frame(a=rbind(A,B,C))
as.bastri(dat.fr)
}

}
