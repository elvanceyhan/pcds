% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PropEdge1D.R
\name{funsMuVarPE1D}
\alias{funsMuVarPE1D}
\alias{mu1PE1D}
\alias{muPE1D}
\alias{fvar1}
\alias{fvar2}
\alias{asyvarPE1D}
\title{Returns the mean and (asymptotic) variance of arc density of Proportional Edge Proximity
Catch Digraph (PE-PCD) for 1D data - middle interval case}
\usage{
mu1PE1D(r, c)

muPE1D(r, c)

fvar1(r, c)

fvar2(r, c)

asyvarPE1D(r, c)
}
\arguments{
\item{r}{A positive real number which serves as the expansion parameter in PE proximity region;
must be \eqn{\ge 1}.}

\item{c}{A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.}
}
\value{
\code{muPE1D} returns the mean and \code{asyvarPE1D} returns the asymptotic variance of the
arc density of PE-PCD for \eqn{U(a,b)} data
}
\description{
The functions \code{muPE1D} and \code{asyvarPE1D} and their auxiliary functions.

\code{muPE1D} returns the mean of the (arc) density of PE-PCD
and \code{asyvarPE1D} returns the (asymptotic) variance of the arc density of PE-PCD
for a given centrality parameter \eqn{c \in (0,1)} and an expansion parameter \eqn{r \ge 1} and for 1D uniform data in a
finite interval \eqn{(a,b)}, i.e., data from \eqn{U(a,b)} distribution.

\code{muPE1D} uses auxiliary (internal) function \code{mu1PE1D} which yields mean (i.e., expected value)
of the arc density of PE-PCD for a given \eqn{c \in (0,1/2)} and \eqn{r \ge 1}.

\code{asyvarPE1D} uses auxiliary (internal) functions \code{fvar1} which yields asymptotic variance
of the arc density of PE-PCD for \eqn{c \in (1/4,1/2)} and \eqn{r \ge 1}; and \code{fvar2} which yields asymptotic variance
of the arc density of PE-PCD for \eqn{c \in (0,1/4)} and \eqn{r \ge 1}.

See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
}
\examples{
\dontrun{
#Examples for muPE1D
muPE1D(1.2,.4)
muPE1D(1.2,.6)

rseq<-seq(1.01,5,by=.1)
cseq<-seq(0.01,.99,by=.1)

lrseq<-length(rseq)
lcseq<-length(cseq)

mu.grid<-matrix(0,nrow=lrseq,ncol=lcseq)
for (i in 1:lrseq)
  for (j in 1:lcseq)
  {
    mu.grid[i,j]<-muPE1D(rseq[i],cseq[j])
  }

persp(rseq,cseq,mu.grid, xlab="r", ylab="c", zlab="mu(r,c)", theta = -30, phi = 30,
expand = 0.5, col = "lightblue", ltheta = 120, shade = 0.05, ticktype = "detailed")
}

\dontrun{
#Examples for asyvarPE1D
asyvarPE1D(1.2,.8)

rseq<-seq(1.01,5,by=.1)
cseq<-seq(0.01,.99,by=.1)

lrseq<-length(rseq)
lcseq<-length(cseq)

var.grid<-matrix(0,nrow=lrseq,ncol=lcseq)
for (i in 1:lrseq)
  for (j in 1:lcseq)
  {
    var.grid[i,j]<-asyvarPE1D(rseq[i],cseq[j])
  }

persp(rseq,cseq,var.grid, xlab="r", ylab="c", zlab="var(r,c)", theta = -30, phi = 30,
expand = 0.5, col = "lightblue", ltheta = 120, shade = 0.05, ticktype = "detailed")
}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{muCS1D}} and \code{\link{asyvarCS1D}}
}
\author{
Elvan Ceyhan
}
