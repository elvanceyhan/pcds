% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassFunctions.R
\name{summary.Extrema}
\alias{summary.Extrema}
\title{Return a summary of a \code{Extrema} \code{object}}
\usage{
\method{summary}{Extrema}(object, ...)
}
\arguments{
\item{object}{An \code{object} of class \code{Extrema}.}

\item{\dots}{Additional parameters for \code{summary}.}
}
\value{
The \code{call} of the \code{object} of class \code{"Extrema"},
the \code{type} of the extrema (i.e. the description
of the extrema), extrema points,
distances from extrema to the reference \code{object}
(e.g. boundary of a triangle),
some of the data points (from which extrema is found).
}
\description{
Returns the below information about the \code{object}:

\code{call} of the function defining the \code{object},
the \code{type} of the extrema (i.e. the description
of the extrema), extrema points,
distances from extrema to the reference \code{object}
(e.g. boundary of a triangle),
some of the data points (from which extrema is found).
}
\examples{
\dontrun{
n<-10
Xp<-runif.std.tri(n)$gen.points
Ext<-cl2edges.std.tri(Xp)
Ext
summary(Ext)
}

}
\seealso{
\code{\link{print.Extrema}},
\code{\link{print.summary.Extrema}},
and \code{\link{plot.Extrema}}
}
