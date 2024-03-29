% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentSim1D.R
\name{num.arcsCSmid.int}
\alias{num.arcsCSmid.int}
\title{Number of Arcs of of Central Similarity Proximity Catch Digraphs (CS-PCDs) - middle interval case}
\usage{
num.arcsCSmid.int(Xp, int, t, c = 0.5)
}
\arguments{
\item{Xp}{A set or \code{vector} of 1D points which constitute the vertices of CS-PCD.}

\item{int}{A \code{vector} of two real numbers representing an interval.}

\item{t}{A positive real number which serves as the expansion parameter in CS proximity region.}

\item{c}{A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
with the default \code{c=.5}.
For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.}
}
\value{
Number of arcs for the PE-PCD whose vertices are the 1D data set, \code{Xp},
with expansion parameter, \eqn{r \ge 1}, and centrality parameter, \eqn{c \in (0,1)}. PE proximity regions are defined only
for \code{Xp} points inside the interval \code{int}, i.e., arcs are possible for such points only.
}
\description{
Returns the number of arcs of of Central Similarity Proximity Catch Digraphs (CS-PCDs) whose vertices are
the given 1D numerical data set, \code{Xp}.

CS proximity region \eqn{N_{CS}(x,t,c)} is defined with respect to the interval
\code{int}\eqn{=(a,b)} for this function. CS proximity region is constructed with expansion parameter \eqn{t>0} and
centrality parameter \eqn{c \in (0,1)}.

Vertex regions are based on the center associated with the centrality parameter \eqn{c \in (0,1)}.
For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)} and for the number of arcs,
loops are not allowed so arcs are only possible for points inside the middle interval \code{int} for this function.

See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
}
\examples{
c<-.4
t<-2
a<-0; b<-10; int<-c(a,b)

n<-10
Xp<-runif(n,a,b)
num.arcsCSmid.int(Xp,int,t,c)

num.arcsCSmid.int(Xp,int,t,c=.3)

num.arcsCSmid.int(Xp,int,t=1.5,c)

#num.arcsCSmid.int(Xp,int,t,c+5) #gives error
#num.arcsCSmid.int(Xp,int,t,c+10)

n<-10  #try also n<-20
Xp<-runif(n,a-5,b+5)
num.arcsCSint(Xp,int,t,c)

Xp<-runif(n,a+10,b+10)
num.arcsCSmid.int(Xp,int,t,c)

n<-10
Xp<-runif(n,a,b)
num.arcsCSmid.int(Xp,int,t,c)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{num.arcsCSend.int}}, \code{\link{num.arcsPEmid.int}}, and \code{\link{num.arcsPEend.int}}
}
\author{
Elvan Ceyhan
}
