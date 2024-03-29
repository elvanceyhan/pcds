% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CentSim2D.R
\name{plotCSarcs.tri}
\alias{plotCSarcs.tri}
\title{The plot of the arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for a
2D data set - one triangle case}
\usage{
plotCSarcs.tri(
  Xp,
  tri,
  t,
  M = c(1, 1, 1),
  asp = NA,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  xlim = NULL,
  ylim = NULL,
  edge.reg = FALSE,
  ...
)
}
\arguments{
\item{Xp}{A set of 2D points which constitute the vertices of the CS-PCD.}

\item{tri}{A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.}

\item{t}{A positive real number which serves as the expansion parameter in CS proximity region.}

\item{M}{A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
which serves as a center in the interior of the triangle \code{tri};
default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.}

\item{asp}{A \code{numeric} value, giving the aspect ratio \eqn{y/x} (default is \code{NA}), see the official help page for \code{asp} by
typing "\code{? asp}".}

\item{main}{An overall title for the plot (default=\code{NULL}).}

\item{xlab, ylab}{Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).}

\item{xlim, ylim}{Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
(default=\code{NULL} for both).}

\item{edge.reg}{A logical argument to add edge regions to the plot, default is \code{edge.reg=FALSE}.}

\item{\dots}{Additional \code{plot} parameters.}
}
\value{
A plot of the arcs of the CS-PCD whose vertices are the points in data set \code{Xp} and the triangle \code{tri}
}
\description{
Plots the arcs of CS-PCD whose vertices are the data points, \code{Xp} and the triangle \code{tri}. CS proximity regions
are constructed with respect to the triangle \code{tri} with expansion parameter \eqn{t>0}, i.e., arcs may exist only
for \code{Xp} points inside the triangle \code{tri}.

Edge regions are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}; default
is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
}
\examples{
\dontrun{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10  #try also n<-20

set.seed(1)
Xp<-runif.tri(n,Tr)$g

M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)

t<-1.5  #try also t<-2

plotCSarcs.tri(Xp,Tr,t,M,main="Arcs of CS-PCD with t=1.5",xlab="",ylab="",edge.reg = TRUE)

# or try the default center
#plotCSarcs.tri(Xp,Tr,t,main="Arcs of CS-PCD with t=1.5",xlab="",ylab="",edge.reg = TRUE);
#M=(arcsCStri(Xp,Tr,r)$param)$c #the part "M=(arcsPEtri(Xp,Tr,r)$param)$cent" is optional,
#for the below annotation of the plot

#can add vertex labels and text to the figure (with edge regions)
txt<-rbind(Tr,M)
xc<-txt[,1]+c(-.02,.02,.02,.03)
yc<-txt[,2]+c(.02,.02,.02,.03)
txt.str<-c("A","B","C","M")
text(xc,yc,txt.str)
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{plotCSarcs}}, \code{\link{plotPEarcs.tri}} and \code{\link{plotASarcs.tri}}
}
\author{
Elvan Ceyhan
}
