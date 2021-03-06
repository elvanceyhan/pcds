% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{rseg.tri}
\alias{rseg.tri}
\title{Generation of points segregated (in a Type I fashion) from the vertices of a triangle}
\usage{
rseg.tri(k, tri, delta)
}
\arguments{
\item{k}{A positive integer representing the number of points to be generated from the segregation pattern
in the triangle, \code{tri}.}

\item{tri}{Three 2D points, stacked row-wise, each row representing a vertex of the triangle.}

\item{delta}{A positive real number in (0,4/9). \code{delta} is the parameter of segregation (that is,
\eqn{\delta 100} \% area around each vertex in each Delaunay triangle is forbidden for point generation).}
}
\value{
A list with the elements
\item{type}{The type of the pattern from which points are to be generated}
\item{mtitle}{The "main" title for the plot of the point pattern}
\item{parameters}{Exclusion parameter, \code{delta}, of the Type I segregation pattern. \code{delta} is in (0,4/9)
\eqn{\delta 100} \% area around each vertex in the triangle \code{tri} is forbidden for point generation.}
\item{ref.points}{The input set of points, i.e., vertices of \code{tri};
reference points, i.e., points from which generated points are segregated.}
\item{gen.points}{The output set of generated points segregated from the vertices of \code{tri}.}
\item{tri.Y}{Logical output, if \code{TRUE} the triangle \code{tri} is also plotted when the
corresponding plot function from the \code{Patterns} object is called.}
\item{desc.pat}{Description of the point pattern}
\item{num.points}{The vector of two numbers, which are the number of generated points
and the number of reference (i.e., vertice of \code{tri}, which is 3 here).}
\item{xlimit,ylimit}{The ranges of the x- and y-coordinates of the reference points, which are the
vertice so the triangle \code{tri}}
}
\description{
An object of class "Patterns".
Generates \code{k} points uniformly in the support for Type I segregation in
a given triangle, \code{tri}.

\code{delta} is the parameter of segregation (that is, \eqn{\delta 100} \% of the area around each vertex
in the triangle is forbidden for point generation). \code{delta} corresponds to \code{eps} in the
standard equilateral triangle \eqn{T_e} as \eqn{delta=4 eps ^2/3} (see \code{rsegTe}  function).

See (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}) for more on
the segregation pattern.
}
\examples{
n<-10  #try also n<-100
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C)
del<-.4

Xdt<-rseg.tri(n,Tr,del)
Xdt
summary(Xdt)
plot(Xdt)

dat<-rseg.tri(n,Tr,del)$g
Xlim<-range(Tr[,1])
Ylim<-range(Tr[,2])
xd<-Xlim[2]-Xlim[1]
yd<-Ylim[2]-Ylim[1]

plot(Tr,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
polygon(Tr)
points(dat)
xc<-Tr[,1]+c(-.01,.01,.01)
yc<-Tr[,2]+c(.02,.02,.02)
txt.str<-c("A","B","C")
text(xc,yc,txt.str)

dat.fr<-data.frame(a=Tr)
rseg.tri(n,dat.fr,del)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{rasc.tri}}, \code{\link{rsegTe}}, \code{\link{rsegIITe}}, and \code{\link{rsegMT}}
}
