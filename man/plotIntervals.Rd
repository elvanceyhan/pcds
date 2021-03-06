% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCDFunctions.r
\name{plotIntervals}
\alias{plotIntervals}
\title{The plot of the subintervals based on \code{Yp} points together with \code{Xp} points}
\usage{
plotIntervals(
  Xp,
  Yp,
  main = "",
  xlab = "",
  ylab = "",
  xlim = NULL,
  ylim = NULL,
  ...
)
}
\arguments{
\item{Xp}{A set of 1D points whose scatter-plot is provided.}

\item{Yp}{A set of 1D points which constitute the end points of the intervals which
partition the real line.}

\item{main}{An overall title for the plot (default="").}

\item{xlab, ylab}{Titles for the x and y axes, respectively (default="" for both).}

\item{xlim, ylim}{Numeric vectors of length 2, giving the x- and y-coordinate ranges
(default=\code{NULL} for both).}

\item{\dots}{Additional \code{plot} parameters.}
}
\value{
Plot of the intervals based on \code{Yp} points and also scatter plot of \code{Xp} points
}
\description{
Plots the \code{Xp} points and the intervals based on \code{Yp} points
points.
}
\examples{
a<-0; b<-10;

#nx is number of X points (target) and ny is number of Y points (nontarget)
nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;

set.seed(1)
Xp<-runif(nx,a,b)
Yp<-runif(ny,a,b)

plotIntervals(Xp,Yp,xlab="",ylab="")
plotIntervals(Xp,Yp+10,xlab="",ylab="")

}
\seealso{
\code{\link{plotPEregsMI}} and \code{\link{plotDeltri}}
}
