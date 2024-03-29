% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentSim2D.R
\name{IarcCSset2pnt.tri}
\alias{IarcCSset2pnt.tri}
\title{The indicator for the presence of an arc from a point in set \code{S} to the point \code{p} for
Central Similarity Proximity Catch Digraphs (CS-PCDs) - one triangle case}
\usage{
IarcCSset2pnt.tri(S, p, tri, t, M = c(1, 1, 1))
}
\arguments{
\item{S}{A set of 2D points. Presence of an arc from a point in \code{S} to point \code{p} is checked
by the function.}

\item{p}{A 2D point. Presence of an arc from a point in \code{S} to point \code{p} is checked
by the function.}

\item{tri}{A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.}

\item{t}{A positive real number which serves as the expansion parameter in CS proximity region
constructed in the triangle \code{tri}.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri};
default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.}
}
\value{
I(\code{p} is in \eqn{\cup_{x in S} N_{CS}(x,t)}), that is, returns 1 if \code{p} is in \code{S} or inside \eqn{N_{CS}(x,t)} for at least
one \eqn{x} in \code{S}, returns 0 otherwise where CS proximity region is constructed with respect to the triangle \code{tri}
}
\description{
Returns I(\code{p} in \eqn{N_{CS}(x,t)} for some \eqn{x} in \code{S}),
that is, returns 1 if \code{p} in \eqn{\cup_{x in S} N_{CS}(x,t)},
returns 0 otherwise.

CS proximity region is constructed with respect to the triangle \code{tri} with
the expansion parameter \eqn{t>0} and edge regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.

Edges of \code{tri}\eqn{=T(A,B,C)}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as edges 3, 1, and 2, respectively.
If \code{p} is not in \code{S} and either \code{p} or all points in \code{S} are outside \code{tri}, it returns 0,
but if \code{p} is in \code{S}, then it always returns 1 regardless of its location (i.e., loops are allowed).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10

set.seed(1)
Xp<-runif.tri(n,Tr)$gen.points

S<-rbind(Xp[1,],Xp[2,])  #try also S<-c(1.5,1)

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)

tau<-.5

IarcCSset2pnt.tri(S,Xp[3,],Tr,tau,M)
IarcCSset2pnt.tri(S,Xp[3,],Tr,t=1,M)
IarcCSset2pnt.tri(S,Xp[3,],Tr,t=1.5,M)

S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
IarcCSset2pnt.tri(S,Xp[3,],Tr,tau,M)
}

}
\seealso{
\code{\link{IarcCSset2pnt.std.tri}}, \code{\link{IarcCStri}}, \code{\link{IarcCSstd.tri}},
\code{\link{IarcASset2pnt.tri}}, and \code{\link{IarcPEset2pnt.tri}}
}
\author{
Elvan Ceyhan
}
