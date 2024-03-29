% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PropEdge2D.R
\name{num.arcsPEtri}
\alias{num.arcsPEtri}
\title{Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs)
and quantities related to the triangle - one triangle case}
\usage{
num.arcsPEtri(Xp, tri, r, M = c(1, 1, 1))
}
\arguments{
\item{Xp}{A set of 2D points
which constitute the vertices of PE-PCD.}

\item{tri}{A \eqn{3 \times 2} matrix with each row
representing a vertex of the triangle.}

\item{r}{A positive real number
which serves as the expansion parameter in PE proximity region;
must be \eqn{\ge 1}.}

\item{M}{A 2D point in Cartesian coordinates
or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri}
or the circumcenter of \code{tri}
which may be entered as "CC" as well;
default is \eqn{M=(1,1,1)}, i.e.,
the center of mass of \code{tri}.}
}
\value{
A \code{list} with the elements
\item{desc}{A short description of the output: number of arcs
and quantities related to the triangle}
\item{num.arcs}{Number of arcs of the PE-PCD}
\item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
\item{ind.in.tri}{The vector of indices of the \code{Xp} points
that reside in the triangle}
\item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
is the support triangle.}
\item{vertices}{Vertices of the digraph, \code{Xp}.}
}
\description{
An object of class \code{"NumArcs"}.
Returns the number of arcs of
Proportional Edge Proximity Catch Digraphs (PE-PCDs)
whose vertices are the
given 2D numerical data set, \code{Xp}.
It also provides number of vertices
(i.e., number of data points inside the triangle)
and indices of the data points that reside in the triangle.

PE proximity region \eqn{N_{PE}(x,r)} is defined
with respect to the triangle, \code{tri}
with expansion parameter \eqn{r \ge 1} and vertex regions are
based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
in the interior of the triangle \code{tri} or
based on circumcenter of \code{tri};
default is \eqn{M=(1,1,1)}, i.e.,
the center of mass of \code{tri}.
For the number of arcs, loops are not
allowed so arcs are only possible for points
inside the triangle \code{tri} for this function.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);

n<-10  #try also n<-20
set.seed(1)
Xp<-runif.tri(n,Tr)$g

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)

Narcs = num.arcsPEtri(Xp,Tr,r=1.25,M)
Narcs
summary(Narcs)
plot(Narcs)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{num.arcsPEstd.tri}}, \code{\link{num.arcsPE}},
\code{\link{num.arcsCStri}}, and \code{\link{num.arcsAStri}}
}
\author{
Elvan Ceyhan
}
