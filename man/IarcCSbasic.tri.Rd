% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentSim2D.R
\name{IarcCSbasic.tri}
\alias{IarcCSbasic.tri}
\title{The indicator for the presence of an arc from a point to another for Central Similarity Proximity Catch
Digraphs (CS-PCDs) - standard basic triangle case}
\usage{
IarcCSbasic.tri(p1, p2, t, c1, c2, M = c(1, 1, 1), re = NULL)
}
\arguments{
\item{p1}{A 2D point whose CS proximity region is constructed.}

\item{p2}{A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
\code{p1} or not.}

\item{t}{A positive real number which serves as the expansion parameter in CS proximity region; must be \eqn{\ge 1}}

\item{c1, c2}{Positive real numbers which constitute the vertex of the standard basic triangle
adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the standard basic triangle or circumcenter of \eqn{T_b};
default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.}

\item{re}{The index of the edge region in \eqn{T_b} containing the point, either \code{1,2,3} or \code{NULL}
(default is \code{NULL}).}
}
\value{
\eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2},
that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
returns 0 otherwise
}
\description{
Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2},
that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{r \ge 1}.

CS proximity region is defined with respect to the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.

Edge regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
barycentric coordinates in the interior of the standard basic triangle \eqn{T_b};
default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
\code{re} is the index of the edge region \code{p1} resides, with default=\code{NULL}.

If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_b}, it returns 0,
but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).

Any given triangle can be mapped to the standard basic triangle
by a combination of rigid body motions (i.e., translation, rotation, and reflection) and scaling,
preserving uniformity of the points in the original triangle. Hence standard basic triangle is useful for simulation
studies under the uniformity hypothesis.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:arc-density-CS;textual}{pcds}).
}
\examples{
\dontrun{
c1<-.4; c2<-.6
A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
Tb<-rbind(A,B,C);

M<-as.numeric(runif.basic.tri(1,c1,c2)$g)

tau<-2

P1<-as.numeric(runif.basic.tri(1,c1,c2)$g)
P2<-as.numeric(runif.basic.tri(1,c1,c2)$g)
IarcCSbasic.tri(P1,P2,tau,c1,c2,M)

P1<-c(.4,.2)
P2<-c(.5,.26)
IarcCSbasic.tri(P1,P2,tau,c1,c2,M)
IarcCSbasic.tri(P1,P1,tau,c1,c2,M)

#or try
Re<-rel.edge.basic.tri(P1,c1,c2,M)$re
IarcCSbasic.tri(P1,P2,tau,c1,c2,M,Re)
IarcCSbasic.tri(P1,P1,tau,c1,c2,M,Re)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{IarcCStri}} and \code{\link{IarcCSstd.tri}}
}
\author{
Elvan Ceyhan
}
