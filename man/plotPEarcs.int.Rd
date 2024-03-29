% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PropEdge1D.R
\name{plotPEarcs.int}
\alias{plotPEarcs.int}
\title{The plot of the arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) for 1D data
(vertices jittered along \eqn{y}-coordinate) - one interval case}
\usage{
plotPEarcs.int(
  Xp,
  int,
  r,
  c = 0.5,
  Jit = 0.1,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  center = FALSE,
  ...
)
}
\arguments{
\item{Xp}{A \code{vector} of 1D points constituting the vertices of the PE-PCD.}

\item{int}{A \code{vector} of two 1D points constituting the end points of the interval.}

\item{r}{A positive real number which serves as the expansion parameter in PE proximity region;
must be \eqn{\ge 1}.}

\item{c}{A positive real number in \eqn{(0,1)} parameterizing the center of the interval
with the default \code{c=.5}.
For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.}

\item{Jit}{A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
\code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis where \code{Jit} equals to
the range of range of \eqn{\{}\code{Xp}, \code{int}\eqn{\}} multiplied by \code{Jit}).}

\item{main}{An overall title for the plot (default=\code{NULL}).}

\item{xlab, ylab}{Titles of the \eqn{x} and \eqn{y} axes in the plot (default=\code{NULL} for both).}

\item{xlim, ylim}{Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
(default=\code{NULL} for both).}

\item{center}{A logical argument, if \code{TRUE}, plot includes the center of the interval \code{int}
as a vertical line in the plot, else center of the interval is not plotted.}

\item{\dots}{Additional \code{plot} parameters.}
}
\value{
A plot of the arcs of PE-PCD whose vertices are the 1D data set \code{Xp} in which vertices are jittered
along \eqn{y}-axis for better visualization.
}
\description{
Plots the arcs of PE-PCD whose vertices are the 1D points, \code{Xp}. PE proximity regions are constructed with
expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)} and the intervals are based on the interval \code{int}\eqn{=(a,b)}
That is, data set \code{Xp}
constitutes the vertices of the digraph and \code{int} determines the end points of the interval.

For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default for \eqn{Jit=.1}) is added to
the \eqn{y}-direction where \code{Jit} equals to the range of \eqn{\{}\code{Xp}, \code{int}\eqn{\}}
multiplied by \code{Jit} with default for \eqn{Jit=.1}).
\code{center} is a logical argument, if \code{TRUE}, plot includes the center of the interval \code{int}
as a vertical line in the plot, else center of the interval is not plotted.

See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
}
\examples{
\dontrun{
r<-2
c<-.4
a<-0; b<-10; int<-c(a,b)

#n is number of X points
n<-10;  #try also n<-20;

set.seed(1)
xf<-(int[2]-int[1])*.1

Xp<-runif(n,a-xf,b+xf)

Xlim=range(Xp,int)
Ylim=.1*c(-1,1)

jit<-.1
set.seed(1)
plotPEarcs.int(Xp,int,r=1.5,c=.3,jit,xlab="",ylab="",center=TRUE)
set.seed(1)
plotPEarcs.int(Xp,int,r=2,c=.3,jit,xlab="",ylab="",center=TRUE)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{plotPEarcs1D}} and \code{\link{plotCSarcs.int}}
}
\author{
Elvan Ceyhan
}
