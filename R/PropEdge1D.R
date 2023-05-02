#PropEdge1D.R;
#Functions for NPE in R^1
#################################################################

#' @title The indicator for the presence of an arc from a point to another for Proportional Edge
#' Proximity Catch Digraphs (PE-PCDs) - middle interval case
#'
#' @description Returns \eqn{I(p_2 \in N_{PE}(p_1,r,c))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{PE}(p_1,r,c)}, returns 0
#' otherwise, where \eqn{N_{PE}(x,r,c)} is the PE proximity region for point \eqn{x} and is constructed with expansion
#' parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)} for the interval \eqn{(a,b)}.
#'
#' PE proximity regions are defined with respect to the middle interval \code{int} and vertex regions are based
#' on the center associated with the centrality parameter \eqn{c \in (0,1)}. For the interval, \code{int}\eqn{=(a,b)}, the
#' parameterized center is \eqn{M_c=a+c(b-a)}. \code{rv} is the index of the vertex region \eqn{p_1} resides, with default=\code{NULL}.
#' If \eqn{p_1} and \eqn{p_2} are distinct and either of them are outside interval \code{int}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations
#' (i.e., loops are allowed in the digraph).
#'
#' See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param p1,x2 1D points; \eqn{p_1} is the point for which the proximity region, \eqn{N_{PE}(p_1,r,c)} is
#' constructed and \eqn{p_2} is the point which the function is checking whether its inside
#' \eqn{N_{PE}(p_1,r,c)} or not.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param rv The index of the vertex region \eqn{p_1} resides, with default=\code{NULL}.
#'
#' @return \eqn{I(p_2 \in N_{PE}(p_1,r,c))} for points \eqn{p_1} and \eqn{p_2} that is, returns 1 if \eqn{p_2} is in \eqn{N_{PE}(p_1,r,c)},
#' returns 0 otherwise
#'
#' @seealso \code{\link{IndNPEend.int}}, \code{\link{IndNCSmid.int}}, and \code{\link{IndNCSend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' IndNPEmid.int(7,5,int,r,c)
#' IndNPEmid.int(1,3,int,r,c)
#'
#' @export IndNPEmid.int
IndNPEmid.int <- function(p1,x2,int,r,c=.5,rv=NULL)
{
  if (!is.point(p1,1) || !is.point(x2,1) )
  {stop('p1 and x2 must be scalars')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (p1==x2 )
  {arc<-1; return(arc); stop}

  y1<-int[1]; y2<-int[2];
  if (p1<y1 || p1>y2 || x2<y1 || x2>y2 )
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  {rv<-rv.mid.int(p1,int,c)$rv #determines the vertex region for 1D point p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2))!=1)
  {stop('vertex index, rv, must be 1 or 2')}}

  arc<-0;
  if (rv==1)
  {
    if ( x2 < y1+r*(p1-y1) ) {arc <-1}
  } else {
    if ( x2 > y2-r*(y2-p1) ) {arc<-1}
  }
  arc
} #end of the function
#'

#################################################################

#' @title Number of Arcs for Proportional Edge Proximity Catch Digraphs (PE-PCDs) - middle interval case
#'
#' @description Returns the number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) whose vertices are the
#' given 1D numerical data set, \code{Xp}. PE proximity region \eqn{N_{PE}(x,r,c)} is defined with respect to the interval
#' \code{int}\eqn{=(a,b)} for this function.
#'
#' PE proximity region is constructed with expansion parameter \eqn{r \ge 1} and
#' centrality parameter \eqn{c \in (0,1)}.
#'
#' Vertex regions are based on the center associated with the centrality parameter \eqn{c \in (0,1)}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)} and for the number of arcs,
#' loops are not allowed so arcs are only possible for points inside the middle interval \code{int} for this function.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return Number of arcs for the PE-PCD whose vertices are the 1D data set, \code{Xp},
#' with expansion parameter, \eqn{r \ge 1}, and centrality parameter, \eqn{c \in (0,1)}. PE proximity regions are defined only
#' for \code{Xp} points inside the interval \code{int}, i.e., arcs are possible for such points only.
#'
#' @seealso \code{\link{NumArcsPEend.int}}, \code{\link{NumArcsPE1D}}, \code{\link{NumArcsCSmid.int}}, and \code{\link{NumArcsCSend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-10
#' Xp<-runif(n,a,b)
#' NumArcsPEmid.int(Xp,int,r,c)
#' NumArcsPEmid.int(Xp,int,r=1.5,c)
#' }
#'
#' @export NumArcsPEmid.int
NumArcsPEmid.int <- function(Xp,int,r,c=.5)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  Xp<-Xp[Xp>=y1 & Xp<=y2] #data points inside the interval int
  n<-length(Xp)
  if (n>0)
  {
    arcs<-0
    for (i in 1:n)
    {x1<-Xp[i]
    v<-rv.mid.int(x1,int,c)$rv
    if (v==1)
    {
      xR<-y1+r*(x1-y1)
      arcs<-arcs+sum(Xp <= xR )-1 #minus 1 is for the loop at dat.int[i]
    } else {
      xL <-y2-r*(y2-x1)
      arcs<-arcs+sum(Xp >= xL)-1 #minus 1 is for the loop at dat.int[i]
    }
    }
  } else
  {arcs<-0}
  arcs
} #end of the function
#'

#################################################################

# funsMuVarPE1D
#'
#' @title Returns the mean and (asymptotic) variance of arc density of Proportional Edge Proximity
#' Catch Digraph (PE-PCD) for 1D data - middle interval case
#'
#' @description
#' The functions \code{muPE1D} and \code{asyvarPE1D} and their auxiliary functions.
#'
#' \code{muPE1D} returns the mean of the (arc) density of PE-PCD
#' and \code{asyvarPE1D} returns the (asymptotic) variance of the arc density of PE-PCD
#' for a given centrality parameter \eqn{c \in (0,1)} and an expansion parameter \eqn{r \ge 1} and for 1D uniform data in a
#' finite interval \eqn{(a,b)}, i.e., data from \eqn{U(a,b)} distribution.
#'
#' \code{muPE1D} uses auxiliary (internal) function \code{mu1PE1D} which yields mean (i.e., expected value)
#' of the arc density of PE-PCD for a given \eqn{c \in (0,1/2)} and \eqn{r \ge 1}.
#'
#' \code{asyvarPE1D} uses auxiliary (internal) functions \code{fvar1} which yields asymptotic variance
#' of the arc density of PE-PCD for \eqn{c \in (1/4,1/2)} and \eqn{r \ge 1}; and \code{fvar2} which yields asymptotic variance
#' of the arc density of PE-PCD for \eqn{c \in (0,1/4)} and \eqn{r \ge 1}.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return \code{muPE1D} returns the mean and \code{asyvarPE1D} returns the asymptotic variance of the
#' arc density of PE-PCD for \eqn{U(a,b)} data
#'
#' @name funsMuVarPE1D
NULL
#'
#' @seealso \code{\link{muCS1D}} and \code{\link{asyvarCS1D}}
#'
#' @rdname funsMuVarPE1D
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
mu1PE1D <- function(r,c)
{
  mean<-0;
  if (r<1/(1-c))
  {
    mean<-r*c^2+1/2*r-r*c;
  } else {
    if (r<1/c)
    {
      mean<-(1/(2*r))*(c^2*r^2+2*r-1-2*c*r);
    } else {
      mean<-(r-1)/r;
    }}
  mean
} #end of the function
#'
#' @rdname funsMuVarPE1D
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for muPE1D
#' muPE1D(1.2,.4)
#' muPE1D(1.2,.6)
#'
#' rseq<-seq(1.01,5,by=.1)
#' cseq<-seq(0.01,.99,by=.1)
#'
#' lrseq<-length(rseq)
#' lcseq<-length(cseq)
#'
#' mu.grid<-matrix(0,nrow=lrseq,ncol=lcseq)
#' for (i in 1:lrseq)
#'   for (j in 1:lcseq)
#'   {
#'     mu.grid[i,j]<-muPE1D(rseq[i],cseq[j])
#'   }
#'
#' persp(rseq,cseq,mu.grid, xlab="r", ylab="c", zlab="mu(r,c)", theta = -30, phi = 30,
#' expand = 0.5, col = "lightblue", ltheta = 120, shade = 0.05, ticktype = "detailed")
#' }
#'
#' @export muPE1D
muPE1D <- function(r,c)
{
  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  mean<-0;
  if (c <= 1/2)
  {
    mean<-mu1PE1D(r,c);
  } else {
    mean<-mu1PE1D(r,1-c);
  }
  mean
} #end of the function
#'
#' @rdname funsMuVarPE1D
#'
fvar1 <- function(r,c)
{
  asyvar<-0;
  if (r<1/(1-c))
  {
    asyvar<--1/3*(12*c^4*r^4-24*c^3*r^4+3*c^2*r^5+15*c^2*r^4-3*c*r^5-9*c^2*r^3-3*c*r^4+r^5+6*c^2*r^2+9*c*r^3-r^4-6*c*r^2-2*r^3+4*r^2-3*r+1)/r^2;
  } else {
    if (r<1/c)
    {
      asyvar<--1/3*(3*c^4*r^4+c^3*r^5-c^3*r^4-11*c^3*r^3-3*c^2*r^4+6*c^2*r^3+9*c^2*r^2+3*c*r^3-9*c*r^2+3*c*r-r^2+2*r-1)/r^2;
    } else {
      asyvar<-1/3*(2*r-3)/r^2;
    }}
  asyvar
} #end of the function
#'
#' @rdname funsMuVarPE1D
#'
fvar2 <- function(r,c)
{
  asyvar<-0;
  if (r<1/(1-c))
  {
    asyvar<--1/3*(12*c^4*r^4-24*c^3*r^4+3*c^2*r^5+15*c^2*r^4-3*c*r^5-9*c^2*r^3-3*c*r^4+r^5+6*c^2*r^2+9*c*r^3-r^4-6*c*r^2-2*r^3+4*r^2-3*r+1)/r^2;
  } else {
    if (r<(1-sqrt(1-4*c))/(2*c))
    {
      asyvar<--1/3*(3*c^4*r^4+c^3*r^5-c^3*r^4-11*c^3*r^3-3*c^2*r^4+6*c^2*r^3+9*c^2*r^2+3*c*r^3-9*c*r^2+3*c*r-r^2+2*r-1)/r^2;
    } else {
      if (r<(1+sqrt(1-4*c))/(2*c))
      {
        asyvar<--1/3*(3*c^4*r^5-c^3*r^5-11*c^3*r^4+3*c^2*r^4+9*c^2*r^3-3*c*r^3-r^2+2*r-1)/r^3;
      } else {
        if (r<1/c)
        {
          asyvar<--1/3*(3*c^4*r^4+c^3*r^5-c^3*r^4-11*c^3*r^3-3*c^2*r^4+6*c^2*r^3+9*c^2*r^2+3*c*r^3-9*c*r^2+3*c*r-r^2+2*r-1)/r^2;
        } else {
          asyvar<-1/3*(2*r-3)/r^2;
        }}}}
  asyvar
} #end of the function
#'
#' @rdname funsMuVarPE1D
#'
#' @examples
#' \dontrun{
#' #Examples for asyvarPE1D
#' asyvarPE1D(1.2,.8)
#'
#' rseq<-seq(1.01,5,by=.1)
#' cseq<-seq(0.01,.99,by=.1)
#'
#' lrseq<-length(rseq)
#' lcseq<-length(cseq)
#'
#' var.grid<-matrix(0,nrow=lrseq,ncol=lcseq)
#' for (i in 1:lrseq)
#'   for (j in 1:lcseq)
#'   {
#'     var.grid[i,j]<-asyvarPE1D(rseq[i],cseq[j])
#'   }
#'
#' persp(rseq,cseq,var.grid, xlab="r", ylab="c", zlab="var(r,c)", theta = -30, phi = 30,
#' expand = 0.5, col = "lightblue", ltheta = 120, shade = 0.05, ticktype = "detailed")
#' }
#'
#' @export asyvarPE1D
asyvarPE1D <- function(r,c)
{
  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  asyvar<-0;
  if (c<1/4)
  {
    asyvar<-fvar2(r,c);
  } else {
    if (c<1/2)
    {
      asyvar<-fvar1(r,c);
    } else {
      if (c<3/4)
      {
        asyvar<-fvar1(r,1-c);
      } else {
        asyvar<-fvar2(r,1-c);
      }}}
  asyvar
} #end of the function
#'

#################################################################

#' @title A test of uniformity of 1D data in a given interval based on Proportional Edge Proximity Catch Digraph
#' (PE-PCD)
#'
#' @description
#' An object of class \code{"htest"}.
#' This is an \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of uniformity of 1D data
#' in one interval based on the normal approximation of the arc density of the PE-PCD with expansion parameter \eqn{r \ge 1}
#' and centrality parameter \eqn{c \in (0,1)}.
#'
#' The function yields the test statistic, \eqn{p}-value for the
#' corresponding \code{alternative}, the confidence interval, estimate and null value for the parameter of interest
#' (which is the arc density), and method and name of the data set used.
#'
#' The null hypothesis is that data is
#' uniform in a finite interval (i.e., arc density of PE-PCD equals to its expected value under uniform
#' distribution) and \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the end
#' points) or right-sided (i.e., data is accumulated around the mid point or center \eqn{M_c}).
#'
#' See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the arc density of PE-PCD based on
#' the 1D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval for the arc density at the given confidence level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate of the parameter, i.e., arc density}
#' \item{null.value}{Hypothesized value for the parameter, i.e., the null arc density, which is usually the
#' mean arc density under uniform distribution.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{TSArcDensCSint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-100  #try also n<-20, 1000
#' Xp<-runif(n,a,b)
#'
#' TSArcDensPEint(Xp,int,r,c)
#' TSArcDensPEint(Xp,int,r,c,alt="g")
#' TSArcDensPEint(Xp,int,r,c,alt="l")
#' }
#'
#' @export TSArcDensPEint
TSArcDensPEint <- function(Xp,int,r,c=.5,alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  Xp<-Xp[Xp>=y1 & Xp<=y2] #data points inside the interval int
  n<-length(Xp)
  if (n<=1)
  {stop('There are not enough Xp points in the interval, int, to compute arc density!')}

  num.arcs<-NumArcsPEint(Xp,int,r,c)$num.arcs
  arc.dens<-num.arcs/(n*(n-1))
  estimate1<-arc.dens
  mn<-muPE1D(r,c)
  asy.var<-asyvarPE1D(r,c)
  TS<-sqrt(n) *(arc.dens-mn)/sqrt(asy.var)
  method <-c("Large Sample z-Test Based on Arc Density of PE-PCD for Testing Uniformity of 1D Data")
  names(estimate1) <-c("arc density")

  null.dens<-mn
  names(null.dens) <-"(expected) arc density"
  names(TS) <-"standardized arc density (i.e., Z)"

  if (alternative == "less") {
    pval <-pnorm(TS)
    cint <-arc.dens+c(-Inf, qnorm(conf.level))*sqrt(asy.var/n)
  }
  else if (alternative == "greater") {
    pval <-pnorm(TS, lower.tail = FALSE)
    cint <-arc.dens+c(-qnorm(conf.level),Inf)*sqrt(asy.var/n)
  }
  else {
    pval <-2 * pnorm(-abs(TS))
    alpha <-1 - conf.level
    cint <-qnorm(1 - alpha/2)
    cint <-arc.dens+c(-cint, cint)*sqrt(asy.var/n)
  }
  attr(cint, "conf.level") <-conf.level

  rval <-list(
    statistic=TS,
    p.value=pval,
    conf.int = cint,
    estimate = estimate1,
    null.value = null.dens,
    alternative = alternative,
    method = method,
    data.name = dname
  )

  attr(rval, "class") <-"htest"
  return(rval)
} #end for the function
#'

#################################################################

#' @title A test of segregation/association based on arc density of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 1D data
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points in the range (i.e., range) of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster away from \code{Yp} points) and association (where \code{Xp} points cluster around
#' \code{Yp} points) based on the normal approximation of the arc density of the PE-PCD for uniform 1D data.
#'
#' The function yields the test statistic, \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval, estimate and null value for the parameter of interest (which is the arc density),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points in the range of \code{Yp} points, arc density
#' of PE-PCD whose vertices are \code{Xp} points equals to its expected value under the uniform distribution and
#' \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated around the centers of the triangles, or segregation).
#'
#' PE proximity region is constructed with the expansion parameter \eqn{r \ge 1} and centrality parameter \code{c} which yields
#' \eqn{M}-vertex regions. More precisely, for a middle interval \eqn{(y_{(i)},y_{(i+1)})}, the center is
#' \eqn{M=y_{(i)}+c(y_{(i+1)}-y_{(i)})} for the centrality parameter \eqn{c \in (0,1)}.
#'
#' **Caveat:** This test is currently a conditional test, where \code{Xp} points are assumed to be random, while \code{Yp} points are
#' assumed to be fixed (i.e., the test is conditional on \code{Yp} points).
#' Furthermore, the test is a large sample test when \code{Xp} points are substantially larger than \code{Yp} points,
#' say at least 5 times more.
#' This test is more appropriate when supports of \code{Xp} and \code{Yp} have a substantial overlap.
#' Currently, the \code{Xp} points outside the range of \code{Yp} points are handled with a range correction (or
#' end interval correction) factor (see the description below and the function code.)
#' However, in the special case of no \code{Xp} points in the range of \code{Yp} points, arc density is taken to be 1,
#' as this is clearly a case of segregation. Removing the conditioning and extending it to the case of non-concurring supports is
#' an ongoing line of research of the author of the package.
#'
#' \code{end.int.cor} is for end interval correction, (default is "no end interval correction", i.e., \code{end.int.cor=FALSE}),
#' recommended when both \code{Xp} and \code{Yp} have the same interval support.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}) for more on the uniformity test based on the arc
#' density of PE-PCDs.
#'
#' @param Xp A set of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set of 1D points which constitute the end points of the partition intervals.
#' @param support.int Support interval \eqn{(a,b)} with \eqn{a<b}. Uniformity of \code{Xp} points in this interval
#' is tested.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number which serves as the centrality parameter in PE proximity region;
#' must be in \eqn{(0,1)} (default \code{c=.5}).
#' @param end.int.cor A logical argument for end interval correction, default is \code{FALSE},
#' recommended when both \code{Xp} and \code{Yp} have the same interval support.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the arc density
#' PE-PCD whose vertices are the 1D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding \code{alternative}.}
#' \item{conf.int}{Confidence interval for the arc density at the given confidence level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate of the parameter, i.e., arc density}
#' \item{null.value}{Hypothesized value for the parameter, i.e., the null arc density, which is usually the
#' mean arc density under uniform distribution.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{TSArcDensPE}}, \code{\link{TSDomPEBin1D}}, and \code{\link{TSArcDensPEint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int=c(a,b)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' TSArcDensPE1D(Xp,Yp,int,r,c)
#' #try also TSArcDensPE1D(Xp,Yp,int,r,c,alt="l") and TSArcDensPE1D(Xp,Yp,int,r,c,alt="g")
#'
#' TSArcDensPE1D(Xp,Yp,int,r,c,end.int.cor = TRUE)
#' }
#'
#' @export TSArcDensPE1D
TSArcDensPE1D <- function(Xp,Yp,support.int,r,c=.5,end.int.cor=FALSE,
                            alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

  if ((!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp))))
  {stop('Xp and Yp must be 1D vectors of numerical entries.')}

  if (length(Yp)<2)
  {stop('Yp must be of length >2')}

  if (!is.point(support.int) || support.int[2]<=support.int[1])
  {stop('support.int must be an interval as (a,b) with a<b')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  Arcs<-NumArcsPE1D(Xp,Yp,r,c)  #uses the default c=.5, unless specified otherwise
  NinR<-Arcs$num.in.range #number of Xp points in the range of Yp points

  num.arcs.ints = Arcs$int.num.arcs #vector of number of arcs in the partition intervals
  n.int = length(num.arcs.ints)
  num.arcs = sum(num.arcs.ints[-c(1,n.int)]) #this is to remove the number of arcs in the end intervals
  num.dat.ints = Arcs$num.in.intervals[-c(1,n.int)] #vector of numbers of data points in the partition intervals
  Wvec<-Arcs$w
  LW<-Wvec/sum(Wvec)

  dat.int.ind = Arcs$data.interval.indices #indices of partition intervals in which data points reside
  mid.ind = which(dat.int.ind!=1 & dat.int.ind!=n.int) #indices of Xp points in range of Yp points (i.e., in middle intervals)
  dat.int.ind = dat.int.ind[mid.ind] #removing the end interval indices
  dat.mid = Xp[mid.ind] #Xp points in range of Yp points (i.e., in middle intervals)
  part.int.mid =  t(Arcs$partition.intervals)[-c(1,n.int),] #middle partition intervals

  ind.Xp1 = which(num.dat.ints==1) #indices of partition intervals containing only one Xp point
  if (length(ind.Xp1)>0)
  {
    for (i in ind.Xp1)
    {
      Xpi = dat.mid[dat.int.ind==i+1]
      int =  part.int.mid[i,]
      npe = NPEint(Xpi,int,r,c)
      num.arcs = num.arcs+(npe[2]-npe[1])/Wvec[i]
    }
  }

  asy.mean0<-muPE1D(r,c)  #asy mean value for the (r,c) pair
  asy.mean<-asy.mean0*sum(LW^2)

  asy.var0<-asyvarPE1D(r,c)  #asy variance value for the (r,c) pair
  asy.var<-asy.var0*sum(LW^3)+4*asy.mean0^2*(sum(LW^3)-(sum(LW^2))^2)

  n<-length(Xp)  #number of X points
  if (NinR  == 0)
  {warning('There is no Xp point in the range of Yp points to compute arc density,
           but as this is clearly a segregation pattern, arc density is taken to be 1!')
    arc.dens=1
    TS0<-sqrt(n)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat
  } else
  {  arc.dens<-num.arcs/(NinR*(NinR-1))
  TS0<-sqrt(NinR)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat}  #arc density
  }
  estimate1<-arc.dens
  estimate2<-asy.mean

  method <-c("Large Sample z-Test Based on Arc Density of PE-PCD for Testing Uniformity of 1D Data ---")
  if (end.int.cor==F)
  {
    TS<-TS0
    method <-c(method, " without End Interval Correction")
  }
  else
  {
    m<-length(Yp)  #number of Y points
    NoutRange<-n-NinR #number of points outside of the range

    prop.out<-NoutRange/n #observed proportion of points outside range
    exp.prop.out<-2/m  #expected proportion of points outside range of Y points

    TS<-TS0+abs(TS0)*sign(prop.out-exp.prop.out)*(prop.out-exp.prop.out)^2
    method <-c(method, " with End Interval Correction")
  }

  names(estimate1) <-c("arc density")
  null.dens<-asy.mean
  names(null.dens) <-"(expected) arc density"
  names(TS) <-"standardized arc density (i.e., Z)"

  if (alternative == "less") {
    pval <-pnorm(TS)
    cint <-arc.dens+c(-Inf, qnorm(conf.level))*sqrt(asy.var/NinR)
  }
  else if (alternative == "greater") {
    pval <-pnorm(TS, lower.tail = FALSE)
    cint <-arc.dens+c(-qnorm(conf.level),Inf)*sqrt(asy.var/NinR)
  }
  else {
    pval <-2 * pnorm(-abs(TS))
    alpha <-1 - conf.level
    cint <-qnorm(1 - alpha/2)
    cint <-arc.dens+c(-cint, cint)*sqrt(asy.var/NinR)
  }

  attr(cint, "conf.level") <-conf.level

  rval <-list(
    statistic=TS,
    p.value=pval,
    conf.int = cint,
    estimate = estimate1,
    null.value = null.dens,
    alternative = alternative,
    method = method,
    data.name = dname
  )

  attr(rval, "class") <-"htest"
  return(rval)
} #end of the function
#'

#################################################################

#' @title The arcs of Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D data - middle intervals case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of PE-PCD.
#'
#' For this function, PE proximity regions are constructed with respect to the intervals
#' based on \code{Yp} points with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points only inside the intervals.
#' It also provides various descriptions and quantities about the arcs of the PE-PCD
#' such as number of arcs, arc density, etc.
#'
#' Vertex regions are based on center \eqn{M_c} of each middle interval.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the intervals.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, they are expansion and centrality parameters.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, i.e., \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of PE-PCD for 1D data in the middle intervals}
#' \item{E}{Heads (or arrow ends) of the arcs of PE-PCD for 1D data in the middle intervals}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{ArcsPEend.int}}, \code{\link{ArcsPE1D}}, \code{\link{ArcsCSmid.int}},
#'  \code{\link{ArcsCSend.int}} and \code{\link{ArcsCS1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10;
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' Arcs<-ArcsPEmid.int(Xp,Yp,r,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' S<-Arcs$S
#' E<-Arcs$E
#'
#' ArcsPEmid.int(Xp,Yp,r,c)
#' ArcsPEmid.int(Xp,Yp+10,r,c)
#'
#' jit<-.1
#' yjit<-runif(nx,-jit,jit)
#'
#' Xlim<-range(Xp,Yp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),
#' main="arcs of PE-PCD for points (jittered along y-axis)\n in middle intervals ",
#' xlab=" ", ylab=" ", xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit),pch=".")
#' abline(h=0,lty=1)
#' points(Xp, yjit,pch=".",cex=3)
#' abline(v=Yp,lty=2)
#' arrows(S, yjit, E, yjit, length = .05, col= 4)
#' }
#'
#' @export ArcsPEmid.int
ArcsPEmid.int <- function(Xp,Yp,r,c=.5)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)

  if (ny<=1 | nx<=1)
  {
    S<-E<-vector(); nx2<-0
  } else
  {
    Xs<-sort(Xp); Ys<-sort(Yp)  #sorted data points from classes X and Y
    ymin<-Ys[1]; ymax<-Ys[ny];

    int<-rep(0,nx)
    for (i in 1:nx)
      int[i]<-(Xs[i]>ymin & Xs[i] < ymax )  #indices of X points in the middle intervals, i.e., inside min(Yp) and max (Yp)

    Xint<-Xs[int==1] # X points inside  min(Yp) and max (Yp)
    XLe<-Xs[Xs<ymin] # X points in the left end interval of Yp points
    XRe<-Xs[Xs>ymax] # X points in the right end interval of Yp points

    nt<-ny-1 #number of Yp middle intervals
    nx2<-length(Xint)  #number of Xp points inside the middle intervals

    if (nx2==0)
    {S<-E<-NA
    } else
    {
      i.int<-rep(0,nx2)
      for (i in 1:nx2)
        for (j in 1:nt)
        {
          if (Xint[i]>=Ys[j] & Xint[i] < Ys[j+1] )
            i.int[i]<-j #indices of the Yp intervals in which X points reside
        }

      #the arcs of PE-PCDs for parameters r and c
      S<-E<-vector()  #S is for source and E is for end points for the arcs for middle intervals
      for (i in 1:nt)
      {
        Xi<-Xint[i.int==i] #X points in the ith Yp mid interval
        ni<-length(Xi)
        if (ni>1 )
        {
          y1<-Ys[i]; y2<-Ys[i+1]; int<-c(y1,y2)
          for (j in 1:ni)
          {x1 <-Xi[j] ; Xinl<-Xi[-j] #to avoid loops
          v<-rv.mid.int(x1,int,c)$rv
          if (v==1)
          {
            xR<-y1+r*(x1-y1)
            ind.tails<-((Xinl < min(xR,y2)) & (Xinl > y1))
            st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
            S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
          } else {
            xL <-y2-r*(y2-x1)
            ind.tails<-((Xinl < y2) & (Xinl > max(xL,y1)))
            st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
            S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
          }
          }
        }

      }
    }
  }
  if (length(S)==0)
  {S<-E<-NA}

  param<-c(r,c)
  names(param)<-c("expansion parameter","centrality parameter")
  typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D Points in the Middle Intervals with Expansion Parameter r = ",r," and Centrality Parameter c = ",c,sep="")

  main.txt<-paste("Arcs of PE-PCD for Points (jittered\n along y-axis) in Middle Intervals with r = ",round(r,2)," and c = ",round(c,2),sep="")

  nvert<-nx2; nint<-ny-1; narcs<-ifelse(sum(is.na(S))==0,length(S),0);
  arc.dens<-ifelse(nvert>1,narcs/(nvert*(nvert-1)),NA)

  quantities<-c(nvert,ny,nint,narcs,arc.dens)
  names(quantities)<-c("number of vertices", "number of partition points",
                       "number of middle intervals","number of arcs", "arc density")

  res<-list(
    type=typ,
    parameters=param,
    tess.points=Yp, tess.name=yname, #tessellation points
    vertices=Xp, vert.name=xname, #vertices of the digraph
    S=S,
    E=E,
    mtitle=main.txt,
    quant=quantities
  )

  class(res)<-"PCDs"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################
#The case of end intervals
#################################################################


#################################################################

#' @title The indicator for the presence of an arc from a point to another for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) - end interval case
#'
#' @description Returns \eqn{I(p_2 \in N_{PE}(p_1,r))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{PE}(p_1,r)}, returns 0
#' otherwise, where \eqn{N_{PE}(x,r)} is the PE proximity region for point \eqn{x} with expansion parameter \eqn{r \ge 1}
#' for the region outside the interval \eqn{(a,b)}.
#'
#' \code{rv} is the index of the end vertex region \eqn{p_1} resides, with default=\code{NULL},
#' and \code{rv=1} for left end interval and \code{rv=2} for the right end interval.
#' If \eqn{p_1} and \eqn{p_2} are distinct and either of them are inside interval \code{int}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param p1 A 1D point whose PE proximity region is constructed.
#' @param p2 A 1D point. The function determines whether \eqn{p_2} is inside the PE proximity region of
#' \eqn{p_1} or not.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param rv Index of the end interval containing the point, either \code{1,2} or \code{NULL} (default is \code{NULL}).
#'
#' @return \eqn{I(p_2 \in N_{PE}(p_1,r))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{PE}(p_1,r)}
#' (i.e., if there is an arc from \eqn{p_1} to \eqn{p_2}), returns 0 otherwise
#'
#' @seealso \code{\link{IndNPEmid.int}}, \code{\link{IndNCSmid.int}}, and \code{\link{IndNCSend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' a<-0; b<-10; int<-c(a,b)
#' r<-2
#'
#' IndNPEend.int(15,17,int,r)
#' IndNPEend.int(1.5,17,int,r)
#' IndNPEend.int(-15,17,int,r)
#'
#' @export IndNPEend.int
IndNPEend.int <- function(p1,p2,int,r,rv=NULL)
{
  if (!is.point(p1,1) || !is.point(p2,1) )
  {stop('p1 and p2 must be scalars')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (p1==p2 )
  {arc<-1; return(arc); stop}

  if ((p1 >= y1 & p1 <= y2) || (p2 >= y1 & p2 <= y2))
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  {rv<-rv.end.int(p1,int)$rv #determines the vertex for the end interval for 1D point p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2))!=1)
  {stop('vertex index, rv, must be 1 or 2')}}

  arc<-0;
  if (rv==1)
  {
    if ( p2 >= y1-r*(y1-p1) & p2< y1 ) {arc <-1}
  } else
  {
    if ( p2 <= y2+r*(p1-y2) & p2>y2 ) {arc<-1}
  }
  arc
} #end of the function
#'

#################################################################

#' @title Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) - end interval case
#'
#' @description Returns the number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) whose
#' vertices are a 1D numerical data set, \code{Xp}, outside the interval \code{int}\eqn{=(a,b)}.
#'
#' PE proximity region is constructed only with expansion parameter \eqn{r \ge 1} for points outside the interval \eqn{(a,b)}.
#' End vertex regions are based on the end points of the interval,
#' i.e., the corresponding vertex region is an interval as \eqn{(-\infty,a)} or \eqn{(b,\infty)} for the interval \eqn{(a,b)}.
#' For the number of arcs, loops are not allowed, so arcs are only possible for points outside
#' the interval, \code{int}, for this function.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A \code{vector} of 1D points which constitute the vertices of the digraph.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return Number of arcs for the PE-PCD with vertices being 1D data set, \code{Xp},
#' expansion parameter, \eqn{r \ge 1}, for the end intervals.
#'
#' @seealso \code{\link{NumArcsPEmid.int}}, \code{\link{NumArcsPE1D}}, \code{\link{NumArcsCSmid.int}}, and \code{\link{NumArcsCSend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-5
#' XpL<-runif(n,a-5,a)
#' XpR<-runif(n,b,b+5)
#' Xp<-c(XpL,XpR)
#'
#' r<-1.2
#' NumArcsPEend.int(Xp,int,r)
#' NumArcsPEend.int(Xp,int,r=2)
#' }
#'
#' @export NumArcsPEend.int
NumArcsPEend.int <- function(Xp,int,r)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  Xp<-Xp[Xp<y1 | Xp>y2]
  n<-length(Xp)
  if (n<=1)
  {arcs<-0
  } else
  {
    arcs<-0
    for (i in 1:n)
    {p1<-Xp[i]; rv<-rv.end.int(p1,int)$rv
    for (j in ((1:n)[-i]) )
    {p2<-Xp[j]
    arcs<-arcs+IndNPEend.int(p1,p2,int,r,rv)
    }
    }
  }
  arcs
} #end of the function
#'

#################################################################

# funsMuVarPEend.int
#'
#' @title Returns the mean and (asymptotic) variance of arc density of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 1D data - end interval case
#'
#' @description
#' Two functions: \code{muPEend.int} and \code{asyvarPEend.int}.
#'
#' \code{muPEend.int} returns the mean of the arc density of PE-PCD
#' and \code{asyvarPEend.int} returns the asymptotic variance of the arc density of PE-PCD
#' for a given expansion parameter \eqn{r \ge 1} for 1D uniform data in the left and right end intervals
#' for the interval \eqn{(a,b)}.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#'
#' @return \code{muPEend.int} returns the mean and \code{asyvarPEend.int} returns the asymptotic variance of the
#' arc density of PE-PCD for uniform data in end intervals
#'
#' @name funsMuVarPEend.int
NULL
#'
#' @seealso \code{\link{muCSend.int}} and \code{\link{asyvarCSend.int}}
#'
#' @rdname funsMuVarPEend.int
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for muPEend.int
#' muPEend.int(1.2)
#'
#' rseq<-seq(1.01,5,by=.1)
#' lrseq<-length(rseq)
#'
#' mu.end<-vector()
#' for (i in 1:lrseq)
#' {
#'   mu.end<-c(mu.end,muPEend.int(rseq[i]))
#' }
#'
#' plot(rseq, mu.end,type="l",
#' ylab=expression(paste(mu,"(r)")),xlab="r",lty=1,xlim=range(rseq),ylim=c(0,1))
#' }
#'
#' @export muPEend.int
muPEend.int <- function(r)
{
  if (!is.point(r,1) || r<1)
  {stop('the argument must be a scalar greater than 1')}

  1-1/(2*r);
} #end of the function
#'
#' @rdname funsMuVarPEend.int
#'
#' @examples
#' \dontrun{
#' #Examples for asyvarPEend.int
#' asyvarPEend.int(1.2)
#'
#' rseq<-seq(1.01,5,by=.1)
#' lrseq<-length(rseq)
#'
#' var.end<-vector()
#' for (i in 1:lrseq)
#' {
#'   var.end<-c(var.end,asyvarPEend.int(rseq[i]))
#' }
#'
#' par(mar=c(5,5,4,2))
#' plot(rseq, var.end,type="l",
#' xlab="r",ylab=expression(paste(sigma^2,"(r)")),lty=1,xlim=range(rseq))
#' }
#'
#' @export asyvarPEend.int
asyvarPEend.int <- function(r)
{
  if (!is.point(r,1) || r<1)
  {stop('the argument must be a scalar greater than 1')}

  (r-1)^2/(3*r^3);
} #end of the function
#'

#################################################################

#' @title The arcs of Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D data - end interval case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of PE-PCD and related parameters and the quantities of the digraph.  \code{Yp} determines the end points of the end intervals.
#'
#' For this function, PE proximity regions are constructed data points outside the intervals based on
#' \code{Yp} points with expansion parameter \eqn{r \ge 1}. That is, for this function,
#' arcs may exist for points only inside end intervals.
#' It also provides various descriptions and quantities about the arcs of the PE-PCD
#' such as number of arcs, arc density, etc.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the intervals.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, it is the expansion parameter.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization based on \code{Yp}.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitutes the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of PE-PCD for 1D data in the end intervals}
#' \item{E}{Heads (or arrow ends) of the arcs of PE-PCD for 1D data in the end intervals}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals (which is 2 for end intervals), number of arcs, and arc density.}
#'
#' @seealso \code{\link{ArcsPEmid.int}}, \code{\link{ArcsPE1D}} , \code{\link{ArcsCSmid.int}},
#' \code{\link{ArcsCSend.int}}  and \code{\link{ArcsCS1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' a<-0; b<-10; int<-c(a,b);
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.5
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)  #try also Yp<-runif(ny,a,b)+c(-10,10)
#'
#' Arcs<-ArcsPEend.int(Xp,Yp,r)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' S<-Arcs$S
#' E<-Arcs$E
#'
#' jit<-.1
#' yjit<-runif(nx,-jit,jit)
#'
#' Xlim<-range(a,b,Xp,Yp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),pch=".",
#' main="arcs of PE-PCDs for points (jittered along y-axis)\n in end intervals ",
#' xlab=" ", ylab=" ", xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit))
#' abline(h=0,lty=1)
#' points(Xp, yjit,pch=".",cex=3)
#' abline(v=Yp,lty=2)
#' arrows(S, yjit, E, yjit, length = .05, col= 4)
#' }
#'
#' @export ArcsPEend.int
ArcsPEend.int <- function(Xp,Yp,r)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))
  rname <-deparse(substitute(r))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  Xs<-sort(Xp); Ys<-sort(Yp)  #sorted data points
  ymin<-Ys[1]; ymax<-max(Yp)

  XLe<-Xs[Xs<ymin]; XRe<-Xs[Xs>ymax] #X points in the left and right end intervals respectively

  #the arcs of PE-PCDs for parameters r and c
  S<-E<-vector()  #S is for source and E is for end points for the arcs

  #for end intervals
  #left end interval
  nle<-length(XLe)
  if (nle>1 )
  {
    for (j in 1:nle)
    {x1 <-XLe[j];  xLe<-XLe[-j] #to avoid loops
    xL<-ymin-r*(ymin-x1)
    ind.tails<-((xLe < ymin) & (xLe > xL))
    st<-sum(ind.tails)  #sum of tails of the arcs with head XLe[j]
    S<-c(S,rep(x1,st)); E<-c(E,xLe[ind.tails])
    }
  }

  #right end interval
  nre<-length(XRe)
  if (nre>1 )
  {

    for (j in 1:nre)
    {x1 <-XRe[j]; xRe<-XRe[-j]
    xR<-ymax+r*(x1-ymax)
    ind.tails<-((xRe < xR) & xRe > ymax )
    st<-sum(ind.tails)  #sum of tails of the arcs with head XRe[j]
    S<-c(S,rep(x1,st)); E<-c(E,xRe[ind.tails])
    }
  }

  if (length(S)==0)
  {S<-E<-NA}

  param<-r
  names(param)<-"expansion parameter"
  typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D Points in the End Intervals with Expansion Parameter r = ",r,sep="")

  main.txt<-paste("Arcs of PE-PCD for Points (jittered\n along y-axis) in End Intervals with r = ",round(r,2),sep="")

  nvert<-nle+nre; ny<-length(Yp); nint<-2; narcs<-ifelse(sum(is.na(S))==0,length(S),0);
  arc.dens<-ifelse(nvert>1,narcs/(nvert*(nvert-1)),NA)

  quantities<-c(nvert,ny,nint,narcs,arc.dens)
  names(quantities)<-c("number of vertices", "number of partition points",
                       "number of end intervals","number of arcs", "arc density")
  res<-list(
    type=typ,
    parameters=param,
    tess.points=Yp, tess.name=yname, #tessellation points
    vertices=Xp, vert.name=xname, #vertices of the digraph
    S=S,
    E=E,
    mtitle=main.txt,
    quant=quantities
  )

  class(res)<-"PCDs"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The plot of the arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) for 1D data
#' (vertices jittered along \eqn{y}-coordinate) - multiple interval case
#'
#' @description Plots the arcs of PE-PCD whose vertices are the 1D points, \code{Xp}. PE proximity regions are constructed with
#' expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)} and the intervals are based on \code{Yp} points (i.e.
#' the intervalization is based on \code{Yp} points). That is, data set \code{Xp}
#' constitutes the vertices of the digraph and \code{Yp} determines the end points of the intervals.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default for \eqn{Jit=.1}) is added to
#' the \eqn{y}-direction where \code{Jit} equals to the range of \code{Xp} and \code{Yp} multiplied by \code{Jit} with default for \eqn{Jit=.1}).
#' \code{centers} is a logical argument, if \code{TRUE}, plot includes the centers of the intervals
#' as vertical lines in the plot, else centers of the intervals are not plotted.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A \code{vector} of 1D points constituting the vertices of the PE-PCD.
#' @param Yp A \code{vector} of 1D points constituting the end points of the intervals.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis
#' where \code{Jit} equals to the range of the union of \code{Xp} and \code{Yp} points multiplied by \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles of the \eqn{x} and \eqn{y} axes in the plot (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param centers A logical argument, if \code{TRUE}, plot includes the centers of the intervals
#' as vertical lines in the plot, else centers of the intervals are not plotted.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of PE-PCD whose vertices are the 1D data set \code{Xp} in which vertices are jittered
#' along \eqn{y}-axis for better visualization.
#'
#' @seealso \code{\link{plotPEarcs.int}} and \code{\link{plotCSarcs1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' Xlim=range(Xp,Yp)
#' Ylim=.1*c(-1,1)
#'
#' jit<-.1
#'
#' set.seed(1)
#' plotPEarcs1D(Xp,Yp,r=1.5,c=.3,jit,xlab="",ylab="",centers=TRUE)
#' set.seed(1)
#' plotPEarcs1D(Xp,Yp,r=2,c=.3,jit,xlab="",ylab="",centers=TRUE)
#' }
#'
#' @export plotPEarcs1D
plotPEarcs1D <- function(Xp,Yp,r,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,centers=FALSE, ...)
{
  arcs<-ArcsPE1D(Xp,Yp,r,c)
  S<-arcs$S
  E<-arcs$E
  nx<-length(Xp)
  ny<-length(Yp)

  if (is.null(xlim))
  {xlim<-range(Xp,Yp)}

  ns<-length(S)
  jit<-(xlim[2]-xlim[1])*Jit
  ifelse(ns<=1,yjit<-0,yjit<-runif(ns,-jit,jit))

  if (is.null(ylim))
  {ylim<-2*c(-jit,jit)}

  if (is.null(main))
  {
    main.text=paste("Arcs of PE-PCD with r = ",r," and c = ",c,sep="")
    if (!centers){
      ifelse(ns<=1,main<-main.text,main<-c(main.text,"\n (arcs jittered along y-axis)"))
    } else
    {
      if (ny<=2)
      {ifelse(ns<=1,main<-c(main.text,"\n (center added)"),main<-c(main.text,"\n (arcs jittered along y-axis & center added)"))
      } else
      {
        ifelse(ns<=1,main<-c(main.text,"\n (centers added)"),main<-c(main.text,"\n (arcs jittered along y-axis & centers added)"))
      }
    }
  }

  plot(Xp, rep(0,nx),main=main, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3, ...)

  if (centers==TRUE)
  {cents<-centersMc(Yp,c)
  abline(v=cents,lty=3,col="green")}
  abline(v=Yp,lty=2,col="red")
  abline(h=0,lty=2)
  if (!is.null(S)) {arrows(S, yjit, E, yjit, length = .05, col= 4)}
} #end of the function
#'

#################################################################
#NPE Functions that work for both middle and end intervals
#################################################################

#' @title The end points of the Proportional Edge (PE) Proximity Region for a point - one interval case
#'
#' @description Returns the end points of the interval which constitutes the PE proximity region for a point in the
#' interval \code{int}\eqn{=(a,b)=}\code{(rv=1,rv=2)}. PE proximity region is constructed with respect to the interval \code{int}
#' with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}.
#'
#' Vertex regions are based on the (parameterized) center, \eqn{M_c},
#' which is \eqn{M_c=a+c(b-a)} for the interval, \code{int}\eqn{=(a,b)}.
#' The PE proximity region is constructed whether \code{x} is inside or outside the interval \code{int}.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param x A 1D point for which PE proximity region is constructed.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return The interval which constitutes the PE proximity region for the point \code{x}
#'
#' @seealso \code{\link{NCSint}}, \code{\link{NPEtri}} and \code{\link{NPEtetra}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' NPEint(7,int,r,c)
#' NPEint(17,int,r,c)
#' NPEint(1,int,r,c)
#' NPEint(-1,int,r,c)
#'
#' @export NPEint
NPEint <- function(x,int,r,c=.5)
{
  if (!is.point(x,1) )
  {stop('x must be a scalar')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>y2)
  {stop('interval is degenerate or void, left end must be smaller than or equal to right end')}

  if (x<y1 || x>y2)
  {
    ifelse(x<y1,reg<-c(y1-r*(y1-x), y1),reg<-c(y2,y2+r*(x-y2)))
  } else
  {
    Mc<-y1+c*(y2-y1)
    if (x<=Mc)
    {
      reg <-c(y1,min(y1+r*(x-y1),y2) )
    } else
    {
      reg<-c(max(y1,y2-r*(y2-x)),y2 )
    }
  }
  reg #proximity region interval
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) - one interval case
#'
#' @description Returns \eqn{I(p_2 \in N_{PE}(p_1,r,c))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{PE}(p_1,r,c)},
#' returns 0 otherwise, where \eqn{N_{PE}(x,r,c)} is the PE proximity region for point \eqn{x} with expansion parameter \eqn{r \ge 1}
#' and centrality parameter \eqn{c \in (0,1)}.
#'
#' PE proximity region is constructed with respect to the
#' interval \eqn{(a,b)}. This function works whether \eqn{p_1} and \eqn{p_2} are inside or outside the interval \code{int}.
#'
#' Vertex regions for middle intervals are based on the center associated with the centrality parameter \eqn{c \in (0,1)}.
#' If \eqn{p_1} and \eqn{p_2} are identical, then it returns 1 regardless of their locations
#' (i.e., loops are allowed in the digraph).
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param p1 A 1D point for which the proximity region is constructed.
#' @param p2 A 1D point for which it is checked whether it resides in the proximity region
#' of \eqn{p_1} or not.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region
#'  must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return \eqn{I(p_2 \in N_{PE}(p_1,r,c))}, that is, returns 1 if \eqn{p_2} in \eqn{N_{PE}(p_1,r,c)}, returns 0 otherwise
#'
#' @seealso \code{\link{IndNPEmid.int}}, \code{\link{IndNPEend.int}} and \code{\link{IndNCSint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' IndNPEint(7,5,int,r,c)
#' IndNPEint(15,17,int,r,c)
#' IndNPEint(1,3,int,r,c)
#'
#' @export IndNPEint
IndNPEint <- function(p1,p2,int,r,c=.5)
{
  if (!is.point(p2,1) )
  {stop('p2 must be a scalar')}

  arc<-0
  pr<-NPEint(p1,int,r,c)  #proximity region as interval
  if (p2>=pr[1] && p2<=pr[2])
  {arc<-1}
  arc
} #end of the function
#'

#################################################################

#' @title The plot of the Proportional Edge (PE) Proximity Regions for a general interval
#' (vertices jittered along \eqn{y}-coordinate) - one interval case
#'
#' @description Plots the points in and outside of the interval \code{int} and also the PE proximity regions (which are also intervals).
#' PE proximity regions are constructed with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default is \eqn{Jit=.1}) times range of proximity
#' regions and \code{Xp}) is added to the \eqn{y}-direction.
#' \code{center} is a logical argument, if \code{TRUE}, plot includes the
#' center of the interval as a vertical line in the plot, else center of the interval is not plotted.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set of 1D points for which PE proximity regions are to be constructed.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis
#' where \code{Jit} equals to the range of the union of \code{Xp} and \code{Yp} points multiplied by \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges.
#' @param center A logical argument, if \code{TRUE}, plot includes the center of the interval
#' as a vertical line in the plot, else center of the interval is not plotted.
#' @param \dots	Additional \code{plot} parameters.
#'
#' @return Plot of the PE proximity regions for 1D points in or outside the interval \code{int}
#'
#' @seealso \code{\link{plotPEregs1D}}, \code{\link{plotCSregs.int}}, and \code{\link{plotCSregs.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-10
#' xf<-(int[2]-int[1])*.1
#' Xp<-runif(n,a-xf,b+xf)  #try also Xp<-runif(n,a-5,b+5)
#' plotPEregs.int(Xp,int,r,c,xlab="x",ylab="")
#'
#' plotPEregs.int(7,int,r,c,xlab="x",ylab="")
#' }
#'
#' @export plotPEregs.int
plotPEregs.int <- function(Xp,int,r,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,center=FALSE, ...)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}
  n<-length(Xp)
  pr<-c()
  for (i in 1:n)
  {x1<-Xp[i]
  pr<-rbind(pr,NPEint(x1,int,r,c))
  }

  if (is.null(xlim))
  {xlim<-range(Xp,int,pr)}

  xr<-xlim[2]-xlim[1]
  jit<-xr*Jit
  ifelse(n<=1,yjit<-0,yjit<-runif(n,-jit,jit))

  if (is.null(ylim))
  {ylim<-2*c(-jit,jit)}

  if (is.null(main))
  {
    main.text=paste("PE Proximity Regions with r = ",r," and c = ",c,sep="")
    if (!center){
      ifelse(n<=1,main<-main.text,main<-c(main.text,"\n (regions jittered along y-axis)"))
    } else
    {
      ifelse(n<=1,main<-c(main.text,"\n (center added)"),main<-c(main.text,"\n (regions jittered along y-axis & center added)"))
    }
  }

  plot(Xp, yjit,main=main, xlab=xlab, ylab=ylab,xlim=xlim+.05*xr*c(-1,1),ylim=ylim,pch=".",cex=3, ...)
  if (center==TRUE)
  {cents<-centersMc(int,c)
  abline(v=cents,lty=3,col="green")}
  abline(v=int,lty=2)
  abline(h=0,lty=2)
  for (i in 1:n)
  {
    plotrix::draw.arc(pr[i,1]+xr*.05,yjit[i],xr*.05, deg1=150,deg2 = 210, col = "blue")
    plotrix::draw.arc(pr[i,2]-xr*.05, yjit[i],xr*.05, deg1=-30,deg2 = 30, col = "blue")
    segments(pr[i,1], yjit[i], pr[i,2], yjit[i], col= "blue")
  }
} #end of the function
#'

#################################################################

#' @title Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) - one interval case
#'
#' @description Returns the number of arcs and various other quantities, vectors, and lists for Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data points in \code{Xp} in the one middle interval case.
#'
#' The data points could be inside or outside the interval is \code{int}\eqn{=(a,b)}. PE proximity region is constructed
#' with an expansion parameter \eqn{r \ge 1} and a centrality parameter \eqn{c \in (0,1)}.
#' \code{int} determines the end points of the interval.
#'
#' The PE proximity region is constructed for both points inside and outside the interval, hence
#' the arcs may exist for all points inside or outside the interval.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set of 1D points which constitute the vertices of PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return A \code{list} with the elements
#' \item{num.arcs}{Total number of arcs in all intervals (including the end intervals),
#' i.e., the number of arcs for the entire PE-PCD}
#' \item{num.in.range}{Number of \code{Xp} points in the interval \code{int}}
#' \item{num.in.ints}{The vector of number of \code{Xp} points in the partition intervals (including the end intervals)}
#' \item{int.num.arcs}{The \code{vector} of the number of arcs of the components of the PE-PCD in the
#' partition intervals (including the end intervals)}
#' \item{data.int.ind}{A \code{vector} of indices of partition intervals in which data points reside.
#' Partition intervals are numbered from left to right with 1 being the left end interval.}
#'
#' @seealso \code{\link{NumArcsPEmid.int}}, \code{\link{NumArcsPEend.int}}, and \code{\link{NumArcsCSint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' xf<-(int[2]-int[1])*.1
#'
#' set.seed(123)
#'
#' n<-10
#' Xp<-runif(n,a-xf,b+xf)
#' NumArcsPEint(Xp,int,r,c)
#' NumArcsPEint(Xp,int,r,c=.3)
#' NumArcsPEint(Xp,int,r=1.5,c)
#' }
#'
#' @export NumArcsPEint
NumArcsPEint <- function(Xp,int,r,c=.5)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  nx<-length(Xp)
  y1<-int[1]; y2<-int[2];

  arcs<-0
  ind.in.tri = NULL
  if (nx<=0)
  {
    arcs<-0
  } else
  {
    int.ind = rep(0,nx)
    ind.mid = which(Xp >= y1 & Xp <= y2)
    dat.mid<-Xp[ind.mid] # Xp points inside the int
    XpL= Xp[Xp<y1]; XpR= Xp[Xp>y2]
    int.ind[which(Xp<y1)]=1
    int.ind[ind.mid]=2
    int.ind[which(Xp>y2)]=3

    # number of arcs for the intervals
    narcs.left = NumArcsPEend.int(XpL,int,r)
    narcs.right = NumArcsPEend.int(XpR,int,r)
    narcs.mid = NumArcsPEmid.int(dat.mid,int,r,c)
    arcs = c(narcs.left,narcs.mid,narcs.right)

    ni.vec = c(length(XpL),length(dat.mid),length(XpR))

    narcs = sum(arcs)
  }

  NinInt = ni.vec[2]

  res<-list(num.arcs=narcs, #number of arcs for the PE-PCD
            num.in.range=NinInt, #number of Xp points in the interval, int
            num.in.intervals=ni.vec, #vector of numbers of Xp points in the partition intervals
            int.num.arcs=arcs, #vector of number of arcs for the partition intervals
            data.interval.indices=int.ind) #indices of partition intervals in which data points reside, i.e., column number of part.int for each Xp point
  res
} #end of the function
#'

#################################################################

#' @title Number of arcs and related quantities of Proportional Edge Proximity Catch Digraphs (PE-PCDs) - multiple interval case
#'
#' @description Returns the number of arcs and various other quantities, vectors, and lists for Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data points in \code{Xp} in the multiple interval case.
#'
#' For this function, PE proximity regions are constructed data points inside or outside the intervals based
#' on \code{Yp} points with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points in the middle or end intervals.
#'
#' Range (or convex hull) of \code{Yp} (i.e., the interval \eqn{(\min(Yp),\max(Yp))}) is partitioned by the spacings based on
#' \code{Yp} points (i.e., multiple intervals are these partition intervals based on the order statistics of \code{Yp} points
#' whose union constitutes the range of \code{Yp} points). For the number of arcs, loops are not counted.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the partition intervals.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside the middle (partition) intervals
#' with the default \code{c=.5}.
#' For an interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{num.arcs}{Total number of arcs in all intervals (including the end intervals),
#' i.e., the number of arcs for the entire PE-PCD}
#' \item{num.in.range}{Number of \code{Xp} points in the range or convex hull of \code{Yp} points}
#' \item{num.in.ints}{The vector of number of \code{Xp} points in the partition intervals (including the end intervals)
#' based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the lengths of the middle partition intervals (i.e., end intervals excluded)
#' based on \code{Yp} points}
#' \item{int.num.arcs}{The \code{vector} of the number of arcs of the components of the PE-PCD in the
#' partition intervals (including the end intervals) based on \code{Yp} points}
#' \item{part.int}{A matrix with columns corresponding to the partition intervals based on \code{Yp} points.}
#' \item{data.int.ind}{A \code{vector} of indices of partition intervals in which data points reside,
#' i.e., column number of \code{part.int} is provided for each \code{Xp} point. Partition intervals are numbered from left to right
#' with 1 being the left end interval.}
#'
#' @seealso \code{\link{NumArcsPEint}}, \code{\link{NumArcsPEmid.int}}, \code{\link{NumArcsPEend.int}}, and \code{\link{NumArcsCS1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b);
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' NumArcsPE1D(Xp,Yp,r,c)
#' NumArcsPE1D(Xp,Yp,r,c=.3)
#' NumArcsPE1D(Xp,Yp,r=1.5,c)
#' }
#'
#' @export NumArcsPE1D
NumArcsPE1D <- function(Xp,Yp,r,c=.5)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)))
  {stop('Xp and Yp must be 1D vector of numerical entries')}

  nx<-length(Xp); ny<-length(Yp)

  if (nx==0 || ny==0)
  {stop('No Xp or Yp points to construct the PE-PCD')}

  Ys<-sort(Yp)  #sorted data points from class  Y
  ymin<-Ys[1]; ymax<-Ys[ny];
  Yrange=c(ymin, ymax)

  int.ind = rep(0,nx)
  dat.mid<-Xp[Xp>=ymin & Xp <= ymax] # Xp points inside  min(Yp) and max (Yp)
  dat.left= Xp[Xp<ymin]; dat.right= Xp[Xp>ymax]
  int.ind[which(Xp<ymin)]=1
  int.ind[which(Xp>ymax)]=ny+1

  #for end intervals
  narcs.left = NumArcsPEend.int(dat.left,Yrange,r)
  narcs.right = NumArcsPEend.int(dat.right,Yrange,r)
  arcs=narcs.left
  narcs<-narcs.left + narcs.right

  #for middle intervals
  n.int<-ny-1 #number of Yp middle intervals
  nx2<-length(dat.mid)  #number of Xp points inside the middle intervals

  Wvec<-Yspacings<-vector()
  for (j in 1:n.int)
  {
    Yspacings = rbind(Yspacings,c(Ys[j],Ys[j+1]))
    Wvec<-c(Wvec,Ys[j+1]-Ys[j])
  }

  part.ints = rbind(c(-Inf,ymin),Yspacings,c(ymax,Inf))

  for (i in 1:n.int) #to determine which interval data points reside
  {
    ind = which(Xp>=Ys[i] & Xp <= Ys[i+1])
    dat.int<-Xp[ind] #X points in the ith Yp mid interval
    int.ind[ind] = i+1
  }

  ni.vec = vector()
  for (i in 1:n.int)
  {
    dat.int<-Xp[int.ind==i+1] #X points in the ith Yp mid interval
    ni.vec = c(ni.vec,length(dat.int))
    narcs.mid = NumArcsPEmid.int(dat.int,Yspacings[i,],r,c)
    arcs = c(arcs,narcs.mid)
    narcs<-narcs + narcs.mid
  }
  ni.vec = c(length(dat.left),ni.vec,length(dat.right))
  arcs = c(arcs,narcs.right) #reordering the number of arcs vector according to the order of the intervals, from left to right

  res<-list(num.arcs=narcs, #number of arcs for the entire PCD
            num.in.range=nx2, #number of Xp points in the range of Yp points
            num.in.intervals=ni.vec, #number of Xp points in the partition intervals
            weight.vec=Wvec, #lengths of the middle partition intervals
            int.num.arcs=arcs, #vector of number of arcs for the partition intervals
            partition.intervals=t(part.ints), #matrix of the partition intervals, each column is one interval
            data.interval.indices=int.ind) #indices of partition intervals in which data points reside, i.e., column number of part.int for each Xp point
  res
} #end of the function
#'

#################################################################

#' @title The arcs of Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D data - one interval case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of PE-PCD.  \code{int} determines the end points of the interval.
#'
#' For this function, PE proximity regions are constructed data points inside or outside the interval based
#' on \code{int} points with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points in the middle or end intervals.
#' It also provides various descriptions and quantities about the arcs of the PE-PCD
#' such as number of arcs, arc density, etc.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the PE-PCD.
#' @param int A \code{vector} of two 1D points which constitutes the end points of the interval.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, they are expansion and centrality parameters.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization of the real line based on \code{int} points.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{int} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of PE-PCD for 1D data}
#' \item{E}{Heads (or arrow ends) of the arcs of PE-PCD for 1D data}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{ArcsPE1D}}, \code{\link{ArcsPEmid.int}}, \code{\link{ArcsPEend.int}}, and \code{\link{ArcsCS1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b);
#'
#' #n is number of X points
#' n<-10; #try also n<-20
#'
#' xf<-(int[2]-int[1])*.1
#'
#' set.seed(1)
#' Xp<-runif(n,a-xf,b+xf)
#'
#' Arcs<-ArcsPEint(Xp,int,r,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#' }
#'
#' @export ArcsPEint
ArcsPEint <- function(Xp,int,r,c=.5)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(int))

  if (!is.point(Xp,length(Xp)) || !is.point(int,length(int)) )
  {stop('Xp and int must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp)
  S<-E<-vector()  #S is for source and E is for end points for the arcs
  if (nx==0)
  {stop('Not enough points to construct PE-PCD')}

  if (nx>1)
  {
    y1<-int[1]; y2<-int[2];

    ind<-rep(0,nx)
    for (i in 1:nx)
      ind[i]<-(Xp[i]>y1 & Xp[i] < y2 )  #indices of X points inside the interval int

    Xint<-Xp[ind==1] # X points inside the interval int
    XLe<-Xp[Xp<y1] # X points in the left end interval of the interval int
    XRe<-Xp[Xp>y2] # X points in the right end interval of the interval int

    #for left end interval
    nle<-length(XLe)
    if (nle>1 )
    {
      for (j in 1:nle)
      {x1 <-XLe[j];  xLe<-XLe[-j] #to avoid loops
      xL<-y1-r*(y1-x1)
      ind.tails<-((xLe < y1) & (xLe > xL))
      st<-sum(ind.tails)  #sum of tails of the arcs with head XLe[j]
      S<-c(S,rep(x1,st)); E<-c(E,xLe[ind.tails])
      }
    }

    #for the middle interval
    nx2<-length(Xint)  #number of Xp points inside the middle interval

    if (nx2>1 )
    {
      for (j in 1:nx2)
      {x1 <-Xint[j] ; Xinl<-Xint[-j] #to avoid loops
      v<-rv.mid.int(x1,int,c)$rv
      if (v==1)
      {
        xR<-y1+r*(x1-y1)
        ind.tails<-((Xinl < min(xR,y2)) & (Xinl > y1))
        st<-sum(ind.tails)  #sum of tails of the arcs with head Xint[j]
        S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
      } else {
        xL <-y2-r*(y2-x1)
        ind.tails<-((Xinl < y2) & (Xinl > max(xL,y1)))
        st<-sum(ind.tails)  #sum of tails of the arcs with head Xint[j]
        S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
      }
      }
    }

    #for right end interval
    nre<-length(XRe)
    if (nre>1 )
    {
      for (j in 1:nre)
      {x1 <-XRe[j]; xRe<-XRe[-j]
      xR<-y2+r*(x1-y2)
      ind.tails<-((xRe < xR) & xRe > y2 )
      st<-sum(ind.tails)  #sum of tails of the arcs with head XRe[j]
      S<-c(S,rep(x1,st)); E<-c(E,xRe[ind.tails])
      }
    }
  }

  if (length(S)==0)
  {S<-E<-NA}

  param<-c(c,r)
  names(param)<-c("centrality parameter","expansion parameter")

  typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D Points with Expansion Parameter r = ",r, " and Centrality Parameter c = ",c,sep="")

  main.txt<-paste("Arcs of PE-PCD with r = ",round(r,2)," and c = ",round(c,2),"\n (arcs jittered along y-axis)",sep="")

  nvert<-nx; nint<-1; narcs<-ifelse(sum(is.na(S))==0,length(S),0);
  arc.dens<-ifelse(nvert>1,narcs/(nvert*(nvert-1)),NA)

  quantities<-c(nvert,2,nint,narcs,arc.dens)
  names(quantities)<-c("number of vertices", "number of partition points",
                       "number of intervals","number of arcs", "arc density")

  res<-list(
    type=typ,
    parameters=param,
    tess.points=int, tess.name=yname, #tessellation points
    vertices=Xp, vert.name=xname, #vertices of the digraph
    S=S,
    E=E,
    mtitle=main.txt,
    quant=quantities
  )

  class(res)<-"PCDs"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The plot of the arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs) for 1D data
#' (vertices jittered along \eqn{y}-coordinate) - one interval case
#'
#' @description Plots the arcs of PE-PCD whose vertices are the 1D points, \code{Xp}. PE proximity regions are constructed with
#' expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)} and the intervals are based on the interval \code{int}\eqn{=(a,b)}
#' That is, data set \code{Xp}
#' constitutes the vertices of the digraph and \code{int} determines the end points of the interval.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default for \eqn{Jit=.1}) is added to
#' the \eqn{y}-direction where \code{Jit} equals to the range of \eqn{\{}\code{Xp}, \code{int}\eqn{\}}
#' multiplied by \code{Jit} with default for \eqn{Jit=.1}).
#' \code{center} is a logical argument, if \code{TRUE}, plot includes the center of the interval \code{int}
#' as a vertical line in the plot, else center of the interval is not plotted.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A \code{vector} of 1D points constituting the vertices of the PE-PCD.
#' @param int A \code{vector} of two 1D points constituting the end points of the interval.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center of the interval
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis where \code{Jit} equals to
#' the range of range of \eqn{\{}\code{Xp}, \code{int}\eqn{\}} multiplied by \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles of the \eqn{x} and \eqn{y} axes in the plot (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param center A logical argument, if \code{TRUE}, plot includes the center of the interval \code{int}
#' as a vertical line in the plot, else center of the interval is not plotted.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of PE-PCD whose vertices are the 1D data set \code{Xp} in which vertices are jittered
#' along \eqn{y}-axis for better visualization.
#'
#' @seealso \code{\link{plotPEarcs1D}} and \code{\link{plotCSarcs.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b)
#'
#' #n is number of X points
#' n<-10;  #try also n<-20;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(n,a-xf,b+xf)
#'
#' Xlim=range(Xp,int)
#' Ylim=.1*c(-1,1)
#'
#' jit<-.1
#' set.seed(1)
#' plotPEarcs.int(Xp,int,r=1.5,c=.3,jit,xlab="",ylab="",center=TRUE)
#' set.seed(1)
#' plotPEarcs.int(Xp,int,r=2,c=.3,jit,xlab="",ylab="",center=TRUE)
#' }
#'
#' @export plotPEarcs.int
plotPEarcs.int <- function(Xp,int,r,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,center=FALSE, ...)
{
  arcs<-ArcsPEint(Xp,int,r,c)
  S<-arcs$S
  E<-arcs$E

  if (is.null(xlim))
  {xlim<-range(Xp,int)}

  jit<-(xlim[2]-xlim[1])*Jit
  ns<-length(S)
  ifelse(ns<=1,{yjit<-0;Lwd=2},{yjit<-runif(ns,-jit,jit);Lwd=1})

  if (is.null(ylim))
  {ylim<-2*c(-jit,jit)}

  if (is.null(main))
  {
    main.text=paste("Arcs of PE-PCD with r = ",r," and c = ",c,sep="")
    if (!center){
      ifelse(ns<=1,main<-main.text,main<-c(main.text,"\n (arcs jittered along y-axis)"))
    } else
    {
      ifelse(ns<=1,main<-c(main.text,"\n (center added)"),main<-c(main.text,"\n (arcs jittered along y-axis & center added)"))
    }
  }

  nx<-length(Xp)
  plot(Xp, rep(0,nx),main=main, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3, ...)
  if (center==TRUE)
  {cents<-centersMc(int,c)
  abline(v=cents,lty=3,col="green")}
  abline(v=int,lty=2,col="red")
  abline(h=0,lty=2)
  arrows(S, yjit, E, yjit, length = .05, col= 4,lwd=Lwd)
} #end of the function
#'

#################################################################

#' @title The arcs of Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D data - multiple interval case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of PE-PCD and related parameters and the quantities of the digraph.
#' \code{Yp} determines the end points of the intervals.
#'
#' For this function, PE proximity regions are constructed data points inside or outside the intervals based
#' on \code{Yp} points with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points in the middle or end intervals.
#' It also provides various descriptions and quantities about the arcs of the PE-PCD
#' such as number of arcs, arc density, etc.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the intervals.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, they are expansion and centrality parameters.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization of the real line based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of PE-PCD for 1D data}
#' \item{E}{Heads (or arrow ends) of the arcs of PE-PCD for 1D data}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{ArcsPEint}}, \code{\link{ArcsPEmid.int}}, \code{\link{ArcsPEend.int}}, and \code{\link{ArcsCS1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b);
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' Arcs<-ArcsPE1D(Xp,Yp,r,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#' }
#'
#' @export ArcsPE1D
ArcsPE1D <- function(Xp,Yp,r,c=.5)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)
  S<-E<-vector()  #S is for source and E is for end points for the arcs
  if (nx==0 || ny==0)
  {stop('Not enough points to construct PE-PCD')}

  if (nx>1)
  {
    Ys<-sort(Yp)  #sorted data points from classes X and Y
    ymin<-Ys[1]; ymax<-Ys[ny];

    int<-rep(0,nx)
    for (i in 1:nx)
      int[i]<-(Xp[i]>ymin & Xp[i] < ymax )  #indices of X points in the middle intervals, i.e., inside min(Yp) and max (Yp)

    Xint<-Xp[int==1] # X points inside  min(Yp) and max (Yp)
    XLe<-Xp[Xp<ymin] # X points in the left end interval of Yp points
    XRe<-Xp[Xp>ymax] # X points in the right end interval of Yp points

    #for left end interval
    nle<-length(XLe)
    if (nle>1 )
    {
      for (j in 1:nle)
      {x1 <-XLe[j];  xLe<-XLe[-j] #to avoid loops
      xL<-ymin-r*(ymin-x1)
      ind.tails<-((xLe < ymin) & (xLe > xL))
      st<-sum(ind.tails)  #sum of tails of the arcs with head XLe[j]
      S<-c(S,rep(x1,st)); E<-c(E,xLe[ind.tails])
      }
    }

    #for middle intervals
    nt<-ny-1 #number of Yp middle intervals
    nx2<-length(Xint)  #number of Xp points inside the middle intervals

    if (nx2>1)
    {
      i.int<-rep(0,nx2)
      for (i in 1:nx2)
        for (j in 1:nt)
        {
          if (Xint[i]>=Ys[j] & Xint[i] < Ys[j+1] )
            i.int[i]<-j #indices of the Yp intervals in which X points reside
        }

      for (i in 1:nt)
      {
        Xi<-Xint[i.int==i] #X points in the ith Yp mid interval
        ni<-length(Xi)
        if (ni>1 )
        {
          y1<-Ys[i]; y2<-Ys[i+1]; int<-c(y1,y2)
          for (j in 1:ni)
          {x1 <-Xi[j] ; Xinl<-Xi[-j] #to avoid loops
          v<-rv.mid.int(x1,int,c)$rv
          if (v==1)
          {
            xR<-y1+r*(x1-y1)
            ind.tails<-((Xinl < min(xR,y2)) & (Xinl > y1))
            st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
            S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
          } else {
            xL <-y2-r*(y2-x1)
            ind.tails<-((Xinl < y2) & (Xinl > max(xL,y1)))
            st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
            S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
          }
          }
        }
      }
    }

    #for right end interval
    nre<-length(XRe)
    if (nre>1 )
    {

      for (j in 1:nre)
      {x1 <-XRe[j]; xRe<-XRe[-j]
      xR<-ymax+r*(x1-ymax)
      ind.tails<-((xRe < xR) & xRe > ymax )
      st<-sum(ind.tails)  #sum of tails of the arcs with head XRe[j]
      S<-c(S,rep(x1,st)); E<-c(E,xRe[ind.tails])
      }
    }
  }
  if (length(S)==0)
  {S<-E<-NA}

  param<-c(c,r)
  names(param)<-c("centrality parameter","expansion parameter")

  typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D Points with Expansion Parameter r = ",r, " and Centrality Parameter c = ",c,sep="")

  main.txt<-paste("Arcs of PE-PCD with r = ",round(r,2)," and c = ",round(c,2),"\n (arcs jittered along y-axis)",sep="")

  nvert<-nx; nint<-ny+1; narcs<-ifelse(sum(is.na(S))==0,length(S),0);
  arc.dens<-ifelse(nvert>1,narcs/(nvert*(nvert-1)),NA)

  quantities<-c(nvert,ny,nint,narcs,arc.dens)
  names(quantities)<-c("number of vertices", "number of partition points",
                       "number of intervals","number of arcs", "arc density")

  res<-list(
    type=typ,
    parameters=param,
    tess.points=Yp, tess.name=yname, #tessellation points
    vertices=Xp, vert.name=xname, #vertices of the digraph
    S=S,
    E=E,
    mtitle=main.txt,
    quant=quantities
  )

  class(res)<-"PCDs"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Incidence matrix for Proportional-Edge Proximity Catch Digraphs (PE-PCDs)
#' for 1D data - multiple interval case
#'
#' @description Returns the incidence matrix for the PE-PCD for a given 1D numerical data set, \code{Xp},
#' as the vertices of the digraph and \code{Yp} determines the end points of the intervals (in the multi-interval case).
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' PE proximity region is constructed
#' with an expansion parameter \eqn{r \ge 1} and a centrality parameter \eqn{c \in (0,1)}.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp a set of 1D points which constitutes the vertices of the digraph.
#' @param Yp a set of 1D points which constitutes the end points of the intervals
#' that partition the real line.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return Incidence matrix for the PE-PCD with vertices being 1D data set, \code{Xp},
#' and \code{Yp} determines the end points of the intervals (in the multi-interval case)
#'
#' @seealso \code{\link{IncMatCS1D}}, \code{\link{IncMatPEtri}}, and \code{\link{IncMatPE}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10;
#' nx<-10; ny<-4
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' IM<-IncMatPE1D(Xp,Yp,r,c)
#' IM
#'
#' dom.greedy(IM)
#' IndUBdom(IM,6)
#' dom.exact(IM)
#' }
#'
#' @export IncMatPE1D
IncMatPE1D <- function(Xp,Yp,r,c=.5)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)
  nt<-ny-1 #number of Yp middle intervals
  if (nx==0 || ny==0)
  {stop('Not enough points to construct PE-PCD')}

  if (nx>=1)
  {
    Ys<-sort(Yp)  #sorted data points from classes X and Y
    ymin<-Ys[1]; ymax<-Ys[ny];

    pr<-c()
    for (i in 1:nx)
    { x1<-Xp[i]
    if (x1<ymin || x1>=ymax)
    {int<-c(ymin,ymax)
    pr<-rbind(pr,NPEint(x1,int,r,c))
    }
    if (nt>=1)
    {
      for (j in 1:nt)
      {
        if (x1>=Ys[j] & x1 < Ys[j+1] )
        {  y1<-Ys[j]; y2<-Ys[j+1]; int<-c(y1,y2)
        pr<-rbind(pr,NPEint(x1,int,r,c))
        }
      }
    }
    }

    inc.mat<-matrix(0, nrow=nx, ncol=nx)
    for (i in 1:nx)
    { reg<-pr[i,]
    for (j in 1:nx)
    {
      inc.mat[i,j]<-sum(Xp[j]>=reg[1] & Xp[j]<=reg[2])
    }
    }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title Incidence matrix for Proportional-Edge Proximity Catch Digraphs (PE-PCDs)
#' for 1D data - one interval case
#'
#' @description Returns the incidence matrix for the PE-PCD for a given 1D numerical data set, \code{Xp},
#' as the vertices of the digraph and \code{int} determines the end points of the interval (in the one interval case).
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' PE proximity region is constructed
#' with an expansion parameter \eqn{r \ge 1} and a centrality parameter \eqn{c \in (0,1)}.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp a set of 1D points which constitutes the vertices of the digraph.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return Incidence matrix for the PE-PCD with vertices being 1D data set, \code{Xp},
#' and \code{int} determines the end points of the intervals (in the one interval case)
#'
#' @seealso \code{\link{IncMatCSint}}, \code{\link{IncMatPE1D}}, \code{\link{IncMatPEtri}}, and \code{\link{IncMatPE}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' xf<-(int[2]-int[1])*.1
#'
#' set.seed(123)
#'
#' n<-10
#' Xp<-runif(n,a-xf,b+xf)
#'
#' IM<-IncMatPEint(Xp,int,r,c)
#' IM
#'
#' dom.greedy(IM)
#' IndUBdom(IM,6)
#' dom.exact(IM)
#'
#' IncMatPEint(Xp,int+10,r,c)
#' }
#'
#' @export IncMatPEint
IncMatPEint <- function(Xp,int,r,c=.5)
{
  if (!is.point(Xp,length(Xp)) )
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); #ny<-length(Yp)

  if (nx==0)
  {stop('Not enough points to construct PE-PCD')}

  if (nx>=1)
  {
   y1=int[1]; y2<-int[2];

    pr<-c() #proximity region
    for (i in 1:nx)
    { x1<-Xp[i]
    pr<-rbind(pr,NPEint(x1,int,r,c))
    }

    inc.mat<-matrix(0, nrow=nx, ncol=nx)
    for (i in 1:nx)
    { reg<-pr[i,]
    for (j in 1:nx)
    {
      inc.mat[i,j]<-sum(Xp[j]>=reg[1] & Xp[j]<=reg[2])
    }
    }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The plot of the Proportional Edge (PE) Proximity Regions (vertices jittered along \eqn{y}-coordinate)
#' - multiple interval case
#'
#' @description Plots the points in and outside of the intervals based on \code{Yp} points and also the PE proximity regions
#' (i.e., intervals). PE proximity region is constructed with expansion parameter \eqn{r \ge 1} and
#' centrality parameter \eqn{c \in (0,1)}.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)}
#' (default is \eqn{Jit=.1}) times range of \code{Xp} and \code{Yp} and the proximity regions (intervals)) is added to the
#' \eqn{y}-direction.
#'
#' \code{centers} is a logical argument, if \code{TRUE},
#' plot includes the centers of the intervals as vertical lines in the plot,
#' else centers of the intervals are not plotted.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set of 1D points for which PE proximity regions are plotted.
#' @param Yp A set of 1D points which constitute the end points of the intervals which
#' partition the real line.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis
#' where \code{Jit} equals to the range of the union of \code{Xp} and \code{Yp} points multiplied by \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param centers A logical argument, if \code{TRUE}, plot includes the centers of the intervals
#' as vertical lines in the plot, else centers of the intervals are not plotted (default is \code{FALSE}).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the PE proximity regions for 1D points located in the middle or end intervals
#' based on \code{Yp} points
#'
#' @seealso \code{\link{plotPEregs1D}}, \code{\link{plotCSregs.int}}, and \code{\link{plotCSregs1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b);
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' plotPEregs1D(Xp,Yp,r,c,xlab="x",ylab="")
#' }
#'
#' @export
plotPEregs1D <- function(Xp,Yp,r,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,centers=FALSE, ...)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)

  if (ny<1 || nx<1)
  {stop('Both Xp and Yp points must be nonempty to construct PE-PCD')}

  LE<-RE<-vector()
  if (nx>=1)
  { Xp<-sort(Xp)
  Ys<-sort(Yp)  #sorted data points from classes X and Y
  ymin<-Ys[1]; ymax<-Ys[ny];

  in.int<-rep(0,nx)
  for (i in 1:nx)
    in.int[i]<-(Xp[i]>ymin & Xp[i] < ymax )  #indices of X points in the middle intervals, i.e., inside min(Yp) and max (Yp)

  Xint<-Xp[in.int==1] # X points inside  min(Yp) and max (Yp)
  XLe<-Xp[Xp<ymin] # X points in the left end interval of Yp points
  XRe<-Xp[Xp>ymax] # X points in the right end interval of Yp points

  #for left end interval
  nle<-length(XLe)
  if (nle>=1 )
  {
    for (j in 1:nle)
    {x1 <-XLe[j]; int<-c(ymin,ymax)
    pr<-NPEint(x1,int,r,c)
    LE<-c(LE,pr[1]); RE<-c(RE,pr[2])
    }
  }

  #for middle intervals
  nt<-ny-1 #number of Yp middle intervals
  nx2<-length(Xint)  #number of Xp points inside the middle intervals

  if (nx2>=1)
  {
    i.int<-rep(0,nx2)
    for (i in 1:nx2)
      for (j in 1:nt)
      {
        if (Xint[i]>=Ys[j] & Xint[i] < Ys[j+1] )
          i.int[i]<-j #indices of the Yp intervals in which X points reside
      }

    for (i in 1:nt)
    {
      Xi<-Xint[i.int==i] #X points in the ith Yp mid interval
      ni<-length(Xi)
      if (ni>=1 )
      {
        y1<-Ys[i]; y2<-Ys[i+1]; int<-c(y1,y2)
        for (j in 1:ni)
        {x1 <-Xi[j] ;
        pr<-NPEint(x1,int,r,c)
        LE<-c(LE,pr[1]); RE<-c(RE,pr[2])
        }
      }

    }
  }

  #for right end interval
  nre<-length(XRe)
  if (nre>=1 )
  {

    for (j in 1:nre)
    {x1 <-XRe[j]; int<-c(ymin,ymax)
    pr<-NPEint(x1,int,r,c)
    LE<-c(LE,pr[1]); RE<-c(RE,pr[2])
    }
  }
  }

  if (is.null(xlim))
  {xlim<-range(Xp,Yp,LE,RE)}

  xr<-xlim[2]-xlim[1]
  jit<-xr*Jit
  ifelse(nx<=1,yjit<-0,yjit<-runif(nx,-jit,jit))

  if (is.null(ylim))
  {ylim<-2*c(-jit,jit)}

  if (is.null(main))
  {
    main.text=paste("PE Proximity Regions with r = ",r," and c = ",c,sep="")
    if (!centers){
      ifelse(nx<=1,main<-main.text,main<-c(main.text,"\n (regions jittered along y-axis)"))
    } else
    {
      if (ny<=2)
      {ifelse(nx<=1,main<-c(main.text,"\n (center added)"),main<-c(main.text,"\n (regions jittered along y-axis & center added)"))
      } else
      {
        ifelse(nx<=1,main<-c(main.text,"\n (centers added)"),main<-c(main.text,"\n (regions jittered along y-axis & centers added)"))
      }
    }
  }

  plot(Xp, yjit,main=main, xlab=xlab, ylab=ylab,xlim=xlim+.05*xr*c(-1,1),ylim=ylim,pch=".",cex=3, ...)
  if (centers==TRUE)
  {cents<-centersMc(Yp,c)
  abline(v=cents,lty=3,col="green")}
  abline(v=Yp,lty=2)
  abline(h=0,lty=2)
  for (i in 1:nx)
  {
    plotrix::draw.arc(LE[i]+xr*.05, yjit[i],xr*.05, deg1=150,deg2 = 210, col = "blue")
    plotrix::draw.arc(RE[i]-xr*.05,  yjit[i],xr*.05, deg1=-30,deg2 = 30, col = "blue")
    segments(LE[i], yjit[i], RE[i], yjit[i], col= "blue")
  }
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for Proportional Edge
#' Proximity Catch Digraphs (PE-PCDs) for an interval
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 1D data set \code{Xp}.
#'
#' PE proximity region is defined with respect to the interval \code{int} with an expansion parameter, \eqn{r \ge 1},
#' and a centrality parameter, \eqn{c \in (0,1)}, so arcs may exist for \code{Xp} points inside the interval \code{int}\eqn{=(a,b)}.
#'
#' Vertex regions are based on the center associated with the centrality parameter \eqn{c \in (0,1)}.
#' \code{rv} is the index of the vertex region \code{p} resides, with default=\code{NULL}.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' @param p A 1D point that is to be tested for being a dominating point or not of the PE-PCD.
#' @param Xp A set of 1D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}; default \code{c=.5}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param rv Index of the vertex region in which the point resides, either \code{1,2} or \code{NULL}
#' (default is \code{NULL}).
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 1D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Gam1PEtri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' r<-2
#' c<-.4
#' a<-0; b<-10
#' int=c(a,b)
#'
#' Mc<-centMc(int,c)
#'
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif(n,a,b)
#'
#' Gam1PEint(Xp[5],Xp,int,r,c)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Gam1PEint(Xp[i],Xp,int,r,c))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' domset<-Xp[ind.gam1]
#' if (length(ind.gam1)==0)
#' {domset<-NA}
#'
#' #or try
#' Rv<-rv.mid.int(Xp[5],int,c)$rv
#' Gam1PEint(Xp[5],Xp,int,r,c,Rv)
#'
#' Xlim<-range(a,b,Xp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),xlab="",pch=".",xlim=Xlim+xd*c(-.05,.05))
#' abline(h=0)
#' points(cbind(Xp,0))
#' abline(v=c(a,b,Mc),col=c(1,1,2),lty=2)
#' points(cbind(domset,0),pch=4,col=2)
#' text(cbind(c(a,b,Mc),-0.1),c("a","b","Mc"))
#'
#' Gam1PEint(2,Xp,int,r,c,ch.data.pnt = FALSE)
#' #gives an error message if ch.data.pnt = TRUE since point p is not a data point in Xp
#' }
#'
#' @export Gam1PEint
Gam1PEint <- function(p,Xp,int,r,c=.5,rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p,1) )
  {stop('p must be a scalar')}

  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,as.matrix(Xp)))
    {stop('p is not a data point in Xp')}
  }

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (is.null(rv))
  {rv<-rv.mid.int(p,int,c)$rv #determines the vertex region for 1D point p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  Xp<-Xp[(Xp>=y1 & Xp<=y2)]
  n<-length(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (IndNPEint(p,Xp[i],int,r,c)==0)
  {dom<-0;}
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

# funsPG2PE1D
#'
#' @title The functions for probability of domination number \eqn{= 2} for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - middle interval case
#'
#' @description
#' The function \code{PG2PE1D} and its auxiliary functions.
#'
#' Returns \eqn{P(\gamma=2)} for PE-PCD whose vertices are a uniform data set of size \code{n} in a finite interval
#' \eqn{(a,b)} where \eqn{\gamma} stands for the domination number.
#'
#' The PE proximity region \eqn{N_{PE}(x,r,c)} is defined with respect to \eqn{(a,b)} with centrality parameter \eqn{c \in (0,1)}
#' and expansion parameter \eqn{r \ge 1}.
#'
#' To compute the probability \eqn{P(\gamma=2)} for PE-PCD in the 1D case,
#' we partition the domain \eqn{(r,c)=(1,\infty) \times (0,1)}, and compute the probability for each partition
#' set. The sample size (i.e., number of vertices or data points) is a positive integer, \code{n}.
#'
#' @section Auxiliary Functions for \code{PG2PE1D}:
#' The auxiliary functions are \code{PG2AI, PG2AII, PG2AIII, PG2AIV, PG2A, PG2Asym, PG2BIII, PG2B, PG2B,
#' PG2Bsym, PG2CIV, PG2C}, and \code{PG2Csym}, each corresponding to a partition of the domain of
#' \code{r} and \code{c}. In particular, the domain partition is handled in 3 cases as
#'
#' CASE A: \eqn{c \in ((3-\sqrt{5})/2, 1/2)}
#'
#' CASE B: \eqn{c \in (1/4,(3-\sqrt{5})/2)} and
#'
#' CASE C: \eqn{c \in (0,1/4)}.
#'
#'
#' @section Case A - \eqn{c \in ((3-\sqrt{5})/2, 1/2)}:
#' In Case A, we compute \eqn{P(\gamma=2)} with
#'
#' \code{PG2AIV(r,c,n)} if \eqn{1 < r < (1-c)/c};
#'
#' \code{PG2AIII(r,c,n)} if \eqn{(1-c)/c< r < 1/(1-c)};
#'
#' \code{PG2AII(r,c,n)} if \eqn{1/(1-c)< r < 1/c};
#'
#' and \code{PG2AI(r,c,n)} otherwise.
#'
#' \code{PG2A(r,c,n)} combines these functions in Case A: \eqn{c \in ((3-\sqrt{5})/2,1/2)}.
#' Due to the symmetry in the PE proximity regions, we use \code{PG2Asym(r,c,n)} for \eqn{c} in
#' \eqn{(1/2,(\sqrt{5}-1)/2)} with the same auxiliary functions
#'
#' \code{PG2AIV(r,1-c,n)} if \eqn{1 < r < c/(1-c)};
#'
#' \code{PG2AIII(r,1-c,n)} if \eqn{(c/(1-c) < r < 1/c};
#'
#' \code{PG2AII(r,1-c,n)} if \eqn{1/c < r < 1/(1-c)};
#'
#' and \code{PG2AI(r,1-c,n)} otherwise.
#'
#' @section Case B - \eqn{c \in (1/4,(3-\sqrt{5})/2)}:
#'
#' In Case B, we compute \eqn{P(\gamma=2)} with
#'
#' \code{PG2AIV(r,c,n)} if \eqn{1 < r < 1/(1-c)};
#'
#' \code{PG2BIII(r,c,n)} if \eqn{1/(1-c) < r < (1-c)/c};
#'
#' \code{PG2AII(r,c,n)} if \eqn{(1-c)/c < r < 1/c};
#'
#' and \code{PG2AI(r,c,n)} otherwise.
#'
#' PG2B(r,c,n) combines these functions in Case B: \eqn{c \in (1/4,(3-\sqrt{5})/2)}.
#' Due to the symmetry in the PE proximity regions, we use \code{PG2Bsym(r,c,n)} for \code{c} in
#' \eqn{((\sqrt{5}-1)/2,3/4)} with the same auxiliary functions
#'
#' \code{PG2AIV(r,1-c,n)} if \eqn{ 1< r < 1/c};
#'
#' \code{PG2BIII(r,1-c,n)} if \eqn{1/c < r < c/(1-c)};
#'
#' \code{PG2AII(r,1-c,n)} if \eqn{c/(1-c) < r < 1/(1-c)};
#'
#' and \code{PG2AI(r,1-c,n)} otherwise.
#'
#' @section Case C - \eqn{c \in (0,1/4)}:
#'
#' In Case C, we compute \eqn{P(\gamma=2)} with
#'
#' \code{PG2AIV(r,c,n)} if \eqn{1< r < 1/(1-c)};
#'
#' \code{PG2BIII(r,c,n)} if \eqn{1/(1-c) < r < (1-\sqrt{1-4 c})/(2 c)};
#'
#' \code{PG2CIV(r,c,n)} if \eqn{(1-\sqrt{1-4 c})/(2 c) < r < (1+\sqrt{1-4 c})/(2 c)};
#'
#' \code{PG2BIII(r,c,n)} if \eqn{(1+\sqrt{1-4 c})/(2 c) < r <1/(1-c)};
#'
#' \code{PG2AII(r,c,n)} if \eqn{1/(1-c) < r < 1/c};
#'
#' and \code{PG2AI(r,c,n)} otherwise.
#'
#' \code{PG2C(r,c,n)} combines these functions in Case C: \eqn{c \in (0,1/4)}.
#' Due to the symmetry in the PE proximity regions, we use \code{PG2Csym(r,c,n)} for \eqn{c \in (3/4,1)}
#' with the same auxiliary functions
#'
#' \code{PG2AIV(r,1-c,n)} if \eqn{1< r < 1/c};
#'
#' \code{PG2BIII(r,1-c,n)} if \eqn{1/c < r < (1-\sqrt{1-4(1-c)})/(2(1-c))};
#'
#' \code{PG2CIV(r,1-c,n)} if \eqn{(1-\sqrt{1-4(1-c)})/(2(1-c)) < r < (1+\sqrt{1-4(1-c)})/(2(1-c))};
#'
#' \code{PG2BIII(r,1-c,n)} if \eqn{(1+\sqrt{1-4(1-c)})/(2(1-c)) < r < c/(1-c)};
#'
#' \code{PG2AII(r,1-c,n)} if \eqn{c/(1-c)< r < 1/(1-c)};
#'
#' and \code{PG2AI(r,1-c,n)} otherwise.
#'
#' Combining Cases A, B, and C, we get our main function \code{PG2PE1D} which computes \eqn{P(\gamma=2)}
#' for any (\code{r,c}) in its domain.
#'
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param n A positive integer representing the size of the uniform data set.
#'
#' @return \eqn{P(}domination number\eqn{\le 1)} for PE-PCD whose vertices are a uniform data set of size \code{n} in a finite
#' interval \eqn{(a,b)}
#'
#' @name funsPG2PE1D
NULL
#'
#' @seealso \code{\link{PG2PEtri}} and  \code{\link{PG2PE1D.asy}}
#'
#' @author Elvan Ceyhan
#'
#' @rdname funsPG2PE1D
#'
PG2AI <- function(r,c,n)
{
  r^2*(2^n*(1/r)^n*r-2^n*(1/r)^n-2*((r-1)/r^2)^n*r)/((r-1)*(r+1)^2);
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2AII <- function(r,c,n)
{
  -1/((r-1)*(r+1)^2)*r*(((r-1)/r^2)^n*r^2-((c*r+1)/r)^n*r^2+(-(c-1)/r)^n*r^2+((c*r^2+c*r-r+1)/r)^n*r+((r-1)*(c*r+c-1)/r)^n+((c*r+1)/r)^n-((c*r^2+c*r-r+1)/r)^n-(-(c-1)/r)^n);
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2AIII <- function(r,c,n)
{
  -1/(r-1)/(r+1)^2*((-(c-1)/r)^n*r^3+(c/r)^n*r^3+(r-1)^n*r^3-(r-1)^(1+n)*r^2+(-(c-1)*r)^n*r^2+(c*r)^n*r^2-(r-1)^n*r^2-r^3+((r-1)*(c*r+c-1)/r)^n*r+(-(r-1)/r*(c*r+c-r))^n*r-r*(-(c-1)/r)^n-r*(c/r)^n-(r-1)^n*r-r^2-(-(c-1)*r)^n-(c*r)^n+(r-1)^n+r+1);
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2AIV <- function(r,c,n)
{
  1/(r-1)/(r+1)^2*(-(-(c-1)/r)^n*r^3-(c/r)^n*r^3-(r-1)^n*r^3+(r-1)^(1+n)*r^2-(-(c-1)*r)^n*r^2-(c*r)^n*r^2+(r-1)^n*r^2+r^3-(-(r-1)/r*(c*r+c-r))^n*r+r*(-(c-1)/r)^n+r*(c/r)^n+(r-1)^n*r+r^2+(-(r-1)*(c*r+c-1))^n+(-(c-1)*r)^n+(c*r)^n-(r-1)^n-r-1);
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2A <- function(r,c,n)
{
  if (r<1)
  {pg2<-0;
  } else {
    if (r<(1-c)/c)
    {
      pg2<-PG2AIV(r,c,n);
    } else {
      if (r<1/(1-c))
      {
        pg2<-PG2AIII(r,c,n);
      } else {
        if (r<1/c)
        {
          pg2<-PG2AII(r,c,n);
        } else {
          pg2<-PG2AI(r,c,n);
        }}}}
  pg2
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2Asym <- function(r,c,n)
{
  if (r<1)
  {pg2<-0;
  } else {
    if (r<c/(1-c))
    {
      pg2<-PG2AIV(r,1-c,n);
    } else {
      if (r<1/c)
      {
        pg2<-PG2AIII(r,1-c,n);
      } else {
        if (r<1/(1-c))
        {
          pg2<-PG2AII(r,1-c,n);
        } else {
          pg2<-PG2AI(r,1-c,n);
        }}}}
  pg2
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2BIII <- function(r,c,n)
{
  -1/(r-1)/(r+1)^2*(r^3*((r-1)/r^2)^n-((c*r+1)/r)^n*r^3+(-(c-1)/r)^n*r^3+((c*r^2+c*r-r+1)/r)^n*r^2+r*((c*r+1)/r)^n-((c*r^2+c*r-r+1)/r)^n*r-r*(-(c-1)/r)^n-(-(r-1)*(c*r+c-1))^n);
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2B <- function(r,c,n)
{
  if (r<1)
  { pg2<-0;
  } else {
    if (r<1/(1-c))
    {
      pg2<-PG2AIV(r,c,n);
    } else {
      if (r<(1-c)/c)
      {
        pg2<-PG2BIII(r,c,n);
      } else {
        if (r<1/c)
        {
          pg2<-PG2AII(r,c,n);
        } else {
          pg2<-PG2AI(r,c,n);
        }}}}
  pg2
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2Bsym <- function(r,c,n)
{
  if (r<1)
  {pg2<-0;
  } else {
    if (r<1/c)
    {
      pg2<-PG2AIV(r,1-c,n);
    } else {
      if (r<c/(1-c))
      {
        pg2<-PG2BIII(r,1-c,n);
      } else {
        if (r<1/(1-c))
        {
          pg2<-PG2AII(r,1-c,n);
        } else {
          pg2<-PG2AI(r,1-c,n);
        }}}}
  pg2
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2CIV <- function(r,c,n)
{
  1/(r+1)*(-r*(-(c-1)/r)^n-c^n*r-c^n+r*((c*r+1)/r)^n);
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2C <- function(r,c,n)
{
  if (r<1)
  { pg2<-0;
  } else {
    if (r<1/(1-c))
    {
      pg2<-PG2AIV(r,c,n);
    } else {
      if (r<(1-sqrt(1-4*c))/(2*c))
      {
        pg2<-PG2BIII(r,c,n);
      } else {
        if (r<(1+sqrt(1-4*c))/(2*c))
        {
          pg2<-PG2CIV(r,c,n);
        } else {
          if (r<(1-c)/c)
          {
            pg2<-PG2BIII(r,c,n);
          } else {
            if (r<1/c)
            {
              pg2<-PG2AII(r,c,n);
            } else {
              pg2<-PG2AI(r,c,n);
            }}}}}}
  pg2
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
PG2Csym <- function(r,c,n)
{
  if (r<1)
  { pg2<-0;
  } else {
    if (r<1/c)
    {
      pg2<-PG2AIV(r,1-c,n);
    } else {
      if (r<(1-sqrt(1-4*(1-c)))/(2*(1-c)))
      {
        pg2<-PG2BIII(r,1-c,n);
      } else {
        if (r<(1+sqrt(1-4*(1-c)))/(2*(1-c)))
        {
          pg2<-PG2CIV(r,1-c,n);
        } else {
          if (r<c/(1-c))
          {
            pg2<-PG2BIII(r,1-c,n);
          } else {
            if (r<1/(1-c))
            {
              pg2<-PG2AII(r,1-c,n);
            } else {
              pg2<-PG2AI(r,1-c,n);
            }}}}}}
  pg2
} #end of the function
#'
#' @rdname funsPG2PE1D
#'
#' @examples
#' #Examples for the main function PG2PE1D
#' r<-2
#' c<-.5
#'
#' PG2PE1D(r,c,n=10)
#' PG2PE1D(r=1.5,c=1/1.5,n=100)
#'
#' @export
PG2PE1D <- function(r,c,n)
{
  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (c<=0 | c>=1)
  { pg2<-0;
  } else {
    if (c < 1/4)
    {
      pg2<-PG2C(r,c,n);
    } else {
      if (c < (3-sqrt(5))/2)
      {
        pg2<-PG2B(r,c,n);
      } else {
        if (c < 1/2)
        {
          pg2<-PG2A(r,c,n);
        } else {
          if (c < (sqrt(5)-1)/2)
          {
            pg2<-PG2Asym(r,c,n);
          } else {
            if (c < 3/4)
            {
              pg2<-PG2Bsym(r,c,n);
            } else {
              pg2<-PG2Csym(r,c,n);
            }}}}}}
  pg2
} #end of the function
#'

#################################################################
#'
#' @title The asymptotic probability of domination number \eqn{= 2} for Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' - middle interval case
#'
#' @description Returns the asymptotic \eqn{P(}domination number\eqn{\le 1)} for PE-PCD whose vertices are a uniform
#' data set in a finite interval \eqn{(a,b)}.
#'
#' The PE proximity region \eqn{N_{PE}(x,r,c)} is defined with respect to \eqn{(a,b)} with centrality parameter \code{c}
#' in \eqn{(0,1)} and expansion parameter \eqn{r=1/\max(c,1-c)}.
#'
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
#' For the interval, \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return The asymptotic \eqn{P(}domination number\eqn{\le 1)} for PE-PCD whose vertices are a uniform data set in a finite
#' interval \eqn{(a,b)}
#'
#' @seealso \code{\link{PG2PE1D}} and \code{\link{PG2PEtri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.5
#'
#' PG2PE1D.asy(c)
#'
#' PG2PE1D.asy(c=1/1.5)
#' PG2PE1D(r=1.5,c=1/1.5,n=10)
#' PG2PE1D(r=1.5,c=1/1.5,n=100)
#'
#' @export PG2PE1D.asy
PG2PE1D.asy <- function(c)
{
  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  rstar<-1/max(c,1-c)  #r value for the non-degenerate asymptotic distribution

  if (c<=0 | c>=1)
  { pg2<-0;
  } else {
    if (c != 1/2)
    {
      pg2<-rstar/(1+rstar);
    } else {
      pg2<-4/9
    }}
  pg2
} #end of the function
#'

#################################################################

#' @title Indices of the intervals where the 1D point(s) reside
#'
#' @description Returns the indices of intervals for all the points in 1D data set, \code{Xp}, as a vector.
#'
#' Intervals are based on \code{Yp} and left end interval is labeled as 1, the next interval as 2, and so on.
#'
#' @param Xp A set of 1D points for which the indices of intervals are to be determined.
#' @param Yp A set of 1D points from which intervals are constructed.
#'
#' @return The \code{vector} of indices of the intervals in which points in the 1D data set, \code{Xp}, reside
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10; int<-c(a,b)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xf<-(int[2]-int[1])*.1
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)  #try also Yp<-runif(ny,a+1,b-1)
#'
#' ind<-ind.int.set(Xp,Yp)
#' ind
#'
#' jit<-.1
#' yjit<-runif(nx,-jit,jit)
#'
#' Xlim<-range(a,b,Xp,Yp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0), xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit),pch=".")
#' abline(h=0)
#' points(Xp, yjit,pch=".",cex=3)
#' abline(v=Yp,lty=2)
#' text(Xp,yjit,labels=factor(ind))
#' }
#'
#' @export ind.int.set
ind.int.set <- function(Xp,Yp)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)))
  {stop('Both arguments must be 1D vectors of numerical entries')}

  nt<-length(Xp)
  ny<-length(Yp)
  ind.set<-rep(0,nt)
  Ys<-sort(Yp)
  ind.set[Xp<Ys[1]]<-1; ind.set[Xp>Ys[ny]]<-ny+1;
  for (i in 1:(ny-1))
  {
    ind<-(Xp>=Ys[i] & Xp<=Ys[i+1] )
    ind.set[ind]<-i+1
  }
  ind.set
} #end of the function
#'

#################################################################

#' @title The domination number of Proportional Edge Proximity Catch Digraph (PE-PCD) for 1D data
#'
#' @description Returns the domination number, a minimum dominating set of PE-PCD whose vertices are the 1D data set \code{Xp},
#' and the domination numbers for partition intervals based on \code{Yp}.
#'
#' \code{Yp} determines the end points of the intervals (i.e., partition the real line via intervalization).
#' It also includes the domination numbers in the end intervals, with interval label 1 for the left end interval
#' and $|Yp|+1$ for the right end interval.
#'
#' PE proximity region is constructed with expansion parameter \eqn{r \ge 1} and centrality parameter \eqn{c \in (0,1)}.
#'
#' @param Xp A set of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set of 1D points which constitute the end points of the intervals which
#' partition the real line.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int} (default \code{c=.5}).
#'
#' @return A \code{list} with three elements
#' \item{dom.num}{Domination number of PE-PCD with vertex set \code{Xp} and expansion parameter \eqn{r \ge 1} and
#' centrality parameter \eqn{c \in (0,1)}.}
#' \item{mds}{A minimum dominating set of the PE-PCD.}
#' \item{ind.mds}{The data indices of the minimum dominating set of the PE-PCD whose vertices are \code{Xp} points.}
#' \item{int.dom.nums}{Domination numbers of the PE-PCD components for the partition intervals.}
#'
#' @seealso \code{\link{PEdom.nd}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10
#' c<-.4
#' r<-2
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' PEdom1D(Xp,Yp,r,c)
#'
#' PEdom1D(Xp,Yp,r,c=.25)
#' PEdom1D(Xp,Yp,r=1.25,c)
#' }
#'
#' @export PEdom1D
PEdom1D <- function(Xp,Yp,r,c=.5)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)))
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp)  #number of Xp points
  ny<-length(Yp)  #number of Yp points

  Ys<-sort(Yp)  #sorted Yp points (ends of the partition intervals)
  nint<-ny+1

  if (nint==0)
  {
    gam<-0; mds<-NULL
  } else
  {
    Int.Ind<-ind.int.set(Xp,Ys)  #indices of intervals in which Xp points in the data fall

    #calculation of the domination number
    gam<-rep(0,nint); mds<-mds.ind<-c()

    for (i in 2:(nint-1)) #2:(nint-1) is to remove the end intervals
    {
      ith.int.ind = Int.Ind==i
      Xpi<-Xp[ith.int.ind]  #points in ith partition interval
      ind.parti = which(ith.int.ind) #indices of Xpi points (wrt to original data indices)

      ni<-length(Xpi)  #number of points in ith interval
      if (ni==0)  #Gamma=0 piece
      {
        gam[i]<-0
      } else
      {
        int<-c(Ys[i-1],Ys[i])  #end points of the ith interval
        cl2Mc = cl2Mc.int(Xpi,int,c) #closest Xp points to the center of ith interval
        Clvert <- as.numeric(cl2Mc$ext)
        Clvert.ind <- cl2Mc$ind # indices of these extrema wrt Xpi

        Ext.ind = ind.parti[Clvert.ind] #indices of these extrema wrt to the original data

        #Gamma=1 piece
        cnt<-0; j<-1;
        while (j<=2 & cnt==0)
        {
          if ( !is.na(Clvert[j]) && Gam1PEint(Clvert[j],Xpi,int,r,c,rv=j)==1)
          {gam[i]<-1; cnt<-1; mds<-c(mds,Clvert[j]); mds.ind=c(mds.ind,Ext.ind[j])
          }
          else
          {j<-j+1}
        }

        #Gamma=2 piece
        if (cnt==0)
        {gam[i]<-2; mds<-c(mds,Clvert); mds.ind=c(mds.ind,Ext.ind)}
      }
    }
  }

  indL = which(Xp<Ys[1]); indR = which(Xp>Ys[ny]); #original data indices in the left and right end intervals
  XpL = Xp[indL]; XpR = Xp[indR]; #data points in the left and right end intervals


  if (length(XpL)>0)
  {
    gamL=1;
    mdsL.ind=which(XpL == min(XpL));  #indices of min in the left data set
    mdsL=XpL[mdsL.ind] #mds set in the left end interval
    ind.mdsL=indL[mdsL.ind]; #data index for the mds of XpL
  } else
  {
  gamL=0;
  mdsL.ind = mdsL = ind.mdsL = NULL
  }

  if (length(XpR)>0)
  {
    gamR=1;
    mdsR.ind=which(XpR == max(XpR));  #indices of max in the right data set
    mdsR=XpR[mdsR.ind] #mds set in the right end interval
    ind.mdsR=indR[mdsR.ind]; #data index for the mds of XpR
  } else
  {
  gamR=0;
  mdsR.ind = mdsR = ind.mdsR = NULL
  }

  mds<-c(mdsL,mds,mdsR)  #a minimum dominating set
  mds.ind=c(ind.mdsL,mds.ind,ind.mdsR)

  gam[1]<-gamL; gam[nint]<-gamR; #c(gamL,gam,gamR) #adding the domination numbers in the end intervals
  Gam<-sum(gam)  #domination number for the entire digraph including the end intervals

  res<- list(dom.num = Gam,  #domination number
             mds = mds, #a minimum dominating set
             ind.mds = mds.ind, #indices of a minimum dominating set (wrt to original data)
             int.dom.nums = gam #domination numbers for the partition intervals
  )

  res
} #end of the function
#'

#################################################################

#' @title The domination number of Proportional Edge Proximity Catch Digraph (PE-PCD) with
#' non-degeneracy centers - multiple interval case
#'
#' @description Returns the domination number, a minimum dominating set of PE-PCD whose vertices are the 1D data set \code{Xp},
#' and the domination numbers for partition intervals based on \code{Yp}
#' when PE-PCD is constructed with vertex regions based on non-degeneracy centers.
#'
#' \code{Yp} determines the end points of the intervals (i.e., partition the real line via intervalization).
#'
#' PE proximity regions are defined with respect to the intervals based on \code{Yp} points with
#' expansion parameter \eqn{r \ge 1} and vertex regions in each interval are based on the centrality parameter \code{c}
#' which is one of the 2 values of \code{c} (i.e., \eqn{c \in \{(r-1)/r,1/r\}})
#' that renders the asymptotic distribution of domination number
#' to be non-degenerate for a given value of \code{r} in \eqn{(1,2)} and \code{c} is center of mass for \eqn{r=2}.
#' These values are called non-degeneracy centrality parameters and the corresponding centers are called
#' nondegeneracy centers.
#'
#' @param Xp A set of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set of 1D points which constitute the end points of the intervals which
#' partition the real line.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,2]} here.
#'
#' @return A \code{list} with three elements
#' \item{dom.num}{Domination number of PE-PCD with vertex set \code{Xp} and expansion parameter \eqn{r in (1,2]} and
#' centrality parameter \eqn{c \in \{(r-1)/r,1/r\}}.}
#' \item{mds}{A minimum dominating set of the PE-PCD.}
#' \item{ind.mds}{The data indices of the minimum dominating set of the PE-PCD whose vertices are \code{Xp} points.}
#' \item{int.dom.nums}{Domination numbers of the PE-PCD components for the partition intervals.}
#'
#' @seealso \code{\link{PEdom.nd}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10
#' r<-1.5
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' PEdom1D.nd(Xp,Yp,r)
#' PEdom1D.nd(Xp,Yp,r=1.25)
#' }
#'
#' @export PEdom1D.nd
PEdom1D.nd <- function(Xp,Yp,r)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)))
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(r,1) || r<= 1 || r>2)
  {stop('r must be a scalar in (1,2]')}

  nx<-length(Xp)  #number of Xp points
  ny<-length(Yp)  #number of Yp points

  Ys<-sort(Yp)  #sorted Yp points (ends of the partition intervals)
  nint<-ny+1

  if (nint==0)
  {
    gam<-0; mds<-NULL
  } else
  {
    Int.Ind<-ind.int.set(Xp,Ys)  #indices of intervals in which Xp points in the data fall

    #calculation of the domination number
    gam<-rep(0,nint); mds<-mds.ind<-c()

    for (i in 1:(nint-2)) #1:(nint-2) is to remove the end intervals
    {
      ith.int.ind = Int.Ind==i+1
      Xpi<-Xp[ith.int.ind]  #points in ith partition interval
      ind.parti = which(ith.int.ind) #indices of Xpi points (wrt to original data indices)

      ni<-length(Xpi)  #number of points in ith interval
      if (ni==0)
      {
        gam[i]<-0
      } else
      { r.c<-sample(1:2,1)  #random center selection from c1,c2
      c.nd<-c((r-1)/r,1/r)[r.c]
      int<-c(Ys[i],Ys[i+1])  #end points of the ith interval
      cl2v = cl2Mc.int(Xpi,int,c.nd)
      Clvert <-as.numeric(cl2v$ext)
      Clvert.ind<-cl2v$ind # indices of these extrema wrt Xpi

      Ext.ind =ind.parti[Clvert.ind] #indices of these extrema wrt to the original data

      #Gamma=1 piece
      cnt<-0; j<-1;
      while (j<=2 & cnt==0)
      {
        if ( !is.na(Clvert[j]) && Gam1PEint(Clvert[j],Xpi,int,r,c.nd,rv=j)==1)
        {gam[i]<-1; cnt<-1; mds<-c(mds,Clvert[j]); mds.ind=c(mds.ind,Ext.ind[j])
        }
        else
        {j<-j+1}
      }

      #Gamma=2 piece
      if (cnt==0)
      {gam[i]<-2; mds<-c(mds,Clvert); mds.ind=c(mds.ind,Ext.ind)}
      }
    }
  }

  indL = which(Xp<Ys[1]); indR = which(Xp>Ys[ny]); #original data indices in the left and right end intervals
  XpL = Xp[indL]; XpR = Xp[indR]; #data points in the left and right end intervals


  if (length(XpL)>0)
  {
    gamL=1;
    mdsL.ind=which(XpL == min(XpL));  #indices of min in the left data set
    mdsL=XpL[mdsL.ind] #mds set in the left end interval
    ind.mdsL=indL[mdsL.ind]; #data index for the mds of XpL
  } else
  {
    gamL=0;
    mdsL.ind = mdsL = ind.mdsL=NULL
  }

  if (length(XpR)>0)
  {
    gamR=1;
    mdsR.ind=which(XpR == max(XpR));  #indices of max in the right data set
    mdsR=XpR[mdsR.ind] #mds set in the right end interval
    ind.mdsR=indR[mdsR.ind]; #data index for the mds of XpR
  } else
  {
    gamR=0;
    mdsR.ind = mdsR = ind.mdsR=NULL
  }

  mds<-c(mdsL,mds,mdsR)  #a minimum dominating set
  mds.ind=c(ind.mdsL,mds.ind,ind.mdsR)

  gam<-c(gamL,gam,gamR) #adding the domination number in the end intervals
  Gam<-sum(gam)  #domination number for the entire digraph including the end intervals

  res<- list(dom.num=Gam,  #domination number
             mds=mds, #a minimum dominating set
             ind.mds =mds.ind, #indices of a minimum dominating set (wrt to original data)
             int.dom.nums=gam #domination numbers for the partition intervals
  )

  res
} #end of the function
#'

#################################################################

#' @title A test of uniformity for 1D data based on domination number of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) - Binomial Approximation
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of
#' uniformity of \code{Xp} points in the support interval \eqn{(a,b)}).
#'
#' The support interval \eqn{(a,b)} is partitioned as \code{(b-a)*(0:nint)/nint}
#' where \code{nint=round(sqrt(nx),0)} and \code{nx} is number of \code{Xp} points, and the test is for testing the uniformity of \code{Xp}
#' points in the interval \eqn{(a,b)}.
#'
#' The null hypothesis is uniformity of \code{Xp} points on \eqn{(a,b)}.
#' The alternative is deviation of distribution of \code{Xp} points from uniformity. The test is based on the (asymptotic) binomial
#' distribution of the domination number of PE-PCD for uniform 1D data in the partition intervals based on partition of \eqn{(a,b)}.
#'
#' The function yields the test statistic, \eqn{p}-value for the corresponding
#' alternative, the confidence interval, estimate and null value for the parameter of interest (which is
#' \eqn{Pr(}domination number\eqn{\le 1)}), and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points in the support interval, probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{\le 1)}) equals to its expected value) and
#' \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the end points of the partition
#' intervals of the support) or right-sided (i.e., data is accumulated around the centers of the partition intervals).
#'
#' PE proximity region is constructed with the expansion parameter \eqn{r \ge 1} and centrality parameter \code{c} which yields
#' \eqn{M}-vertex regions. More precisely \eqn{M_c=a+c(b-a)} for the centrality parameter \code{c} and for a given \eqn{c \in (0,1)}, the
#' expansion parameter \eqn{r} is taken to be \eqn{1/\max(c,1-c)} which yields non-degenerate asymptotic distribution of the
#' domination number.
#'
#' The test statistic is based on the binomial distribution, when success is defined as domination number being less than
#' or equal to 1 in the one interval case (i.e., number of failures is equal to number of times restricted domination number = 1
#' in the intervals).
#' That is, the test statistic is based on the domination number for \code{Xp} points inside the partition intervals
#' for the PE-PCD. For this approach to work, \code{Xp} must be large for each partition interval,
#' but 5 or more per partition interval seems to work in practice.
#'
#' Probability of success is chosen in the following way for various parameter choices.
#' \code{asy.bin} is a logical argument for the use of asymptotic probability of success for the binomial distribution,
#' default is \code{asy.bin=FALSE}.  When \code{asy.bin=TRUE}, asymptotic probability of success for the binomial distribution is used.
#' When \code{asy.bin=FALSE}, the finite sample probability of success for the binomial distribution is used with number
#' of trials equals to expected number of \code{Xp} points per partition interval.
#'
#' @param Xp A set of 1D points which constitute the vertices of the PE-PCD.
#' @param support.int Support interval \eqn{(a,b)} with \eqn{a<b}. Uniformity of \code{Xp} points in this interval
#' is tested.
#' @param c A positive real number which serves as the centrality parameter in PE proximity region;
#' must be in \eqn{(0,1)} (default \code{c=.5}).
#' @param asy.bin A logical argument for the use of asymptotic probability of success for the binomial distribution,
#' default \code{asy.bin=FALSE}. When \code{asy.bin=TRUE}, asymptotic probability of success for the binomial distribution is used.
#' When \code{asy.bin=FALSE}, the finite sample asymptotic probability of success for the binomial distribution is used with number
#' of trials equals to expected number of \code{Xp} points per partition interval.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{\le 1)} for PE-PCD whose vertices are the 1D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval for \eqn{Pr(}domination number\eqn{\le 1)} at the given level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{A \code{vector} with two entries: first is is the estimate of the parameter, i.e.,
#' \eqn{Pr(}domination number\eqn{\le 1)} and second is the domination number}
#' \item{null.value}{Hypothesized value for the parameter, i.e., the null value for \eqn{Pr(}domination number\eqn{\le 1)}}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{TSDomPEBin}}, \code{\link{PEdom1D}} and \code{\link{PEdom1D.nd}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10; supp<-c(a,b)
#' c<-.4
#'
#' r<-1/max(c,1-c)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#'
#' TSDomPEBin1Dint(Xp,supp,c,alt="t")
#' TSDomPEBin1Dint(Xp,support.int = supp,c=c,alt="t")
#' TSDomPEBin1Dint(Xp,supp,c,alt="l")
#' TSDomPEBin1Dint(Xp,supp,c,alt="g")
#' TSDomPEBin1Dint(Xp,supp,c,alt="t",asy.bin = TRUE)
#' }
#'
#' @export TSDomPEBin1Dint
TSDomPEBin1Dint <- function(Xp,support.int,c=.5,asy.bin=FALSE,
                            alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(quote(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be 1D vectors of numerical entries.')}

  if (!is.point(support.int) || support.int[2]<=support.int[1])
  {stop('support.int must be an interval as (a,b) with a<b.')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  rstar<-1/max(c,1-c)  #r value for the non-degenerate asymptotic distribution for a given c

  nx<-length(Xp)  #number of Xp points
  nint<-round(sqrt(nx),0)
  Yp<-support.int[1]+(support.int[2]-support.int[1])*(0:nint)/nint #Y points (ends of the partition intervals)

  if (asy.bin==TRUE)
  {p<-PG2PE1D.asy(c)
  method <-c("Binomial Test based on the Domination Number for Testing Uniformity of 1D Data \n
             (using the asymptotic probability of success in the binomial distribution)")
  } else
  {
    Enx<-nx/nint
    p<-PG2PE1D(rstar,c,Enx)  #p: prob of success; on average n/nint X points fall on each interval
    method <-c("Binomial Test based on the Domination Number for Testing Uniformity of 1D Data \n
             (using the finite sample probability of success in the binomial distribution)")
  }

  if (length(Yp)<2)
  {stop('Number of partition intervals must be of length >2')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  ny<-length(Yp)
  nint<-ny-1  #number of partition intervals

  dom.num=PEdom1D(Xp,Yp,rstar,c)
  Gammas<- dom.num$int.dom.nums #vector of domination numbers of the partition intervals
  estimate2<- dom.num$dom
  # Bm<-Gam-nint; #the binomial test statistic
  Bm<-sum(Gammas<=1) #sum((3-Gammas)[ind0]>0) #sum(Gammas-2>0); #the binomial test statistic, success is dom num <= 1

  x<-Bm
  pval <-switch(alternative, less = pbinom(x, nint, p),
                greater = pbinom(x - 1, nint, p, lower.tail = FALSE),
                two.sided = {if (p == 0) (x == 0) else if (p == 1) (x == nint)
                  else { relErr <-1 + 1e-07
                  d <-dbinom(x, nint, p)
                  m <-nint * p
                  if (x == m) 1 else if (x < m)
                  {i <-seq.int(from = ceiling(m), to = nint)
                  y <-sum(dbinom(i, nint, p) <= d * relErr)
                  pbinom(x, nint, p) + pbinom(nint - y, nint, p, lower.tail = FALSE)
                  } else {
                    i <-seq.int(from = 0, to = floor(m))
                    y <-sum(dbinom(i, nint, p) <= d * relErr)
                    pbinom(y - 1, nint, p) + pbinom(x - 1, nint, p, lower.tail = FALSE)
                  }
                  }
                })

  p.L <- function(x, alpha) {
    if (x == 0)
      0
    else qbeta(alpha, x, nint - x + 1)
  }
  p.U <- function(x, alpha) {
    if (x == nint)
      1
    else qbeta(1 - alpha, x + 1, nint - x)
  }
  cint <-switch(alternative, less = c(0, p.U(x, 1 - conf.level)),
                greater = c(p.L(x, 1 - conf.level), 1), two.sided = {
                  alpha <-(1 - conf.level)/2
                  c(p.L(x, alpha), p.U(x, alpha))
                })
  attr(cint, "conf.level") <-conf.level

  estimate1 <-x/nint
  names(x) <-"#(domination number is <= 1)" #"domination number - number of partition intervals"
  names(nint) <-"number of partition intervals based on Yp"
  names(p) <-"Pr(Domination Number <= 1)"
  names(estimate1) <-c(" Pr(domination number <= 1)")
  names(estimate2) <-c("|| domination number")

  structure(
    list(statistic = x,
         p.value = pval,
         conf.int = cint,
         estimate = c(estimate1,estimate2),
         null.value = p,
         alternative = alternative,
         method = method,
         data.name = dname
    ),
    class = "htest")

} #end of the function
#'

#################################################################

#' @title A test of segregation/association based on domination number of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 1D data - Binomial Approximation
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points within the partition intervals based on \code{Yp} points (both residing in the
#' support interval \eqn{(a,b)}).
#' The test is for testing the spatial interaction between \code{Xp} and \code{Yp} points.
#'
#' The null hypothesis is uniformity of \code{Xp} points on \eqn{(y_{\min},y_{\max})} (by default)
#' where \eqn{y_{\min}} and \eqn{y_{\max}} are minimum and maximum of \code{Yp} points, respectively.
#' \code{Yp} determines the end points of the intervals (i.e., partition the real line via its spacings called intervalization)
#' where end points are the order statistics of \code{Yp} points.
#'
#' The alternatives are segregation (where \code{Xp} points cluster away from \code{Yp} points i.e., cluster around the centers of the
#' partition intervals) and association (where \code{Xp} points cluster around \code{Yp} points). The test is based on the (asymptotic) binomial
#' distribution of the domination number of PE-PCD for uniform 1D data in the partition intervals based on \code{Yp} points.
#'
#' The test by default is restricted to the range of \code{Yp} points, and so ignores \code{Xp} points outside this range.
#' However, a correction for the \code{Xp} points outside the range of \code{Yp} points is available by setting
#' \code{end.int.cor=TRUE}, which is recommended when both \code{Xp} and \code{Yp} have the same interval support.
#'
#' The function yields the test statistic, \eqn{p}-value for the corresponding
#' alternative, the confidence interval, estimate and null value for the parameter of interest (which is
#' \eqn{Pr(}domination number\eqn{\le 1)}), and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points in the intervals based on \code{Yp} points, probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{\le 1)}) equals to its expected value) and
#' \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated around the centers of the partition intervals, or segregation).
#'
#' PE proximity region is constructed with the expansion parameter \eqn{r \ge 1} and centrality parameter \code{c} which yields
#' \eqn{M}-vertex regions. More precisely, for a middle interval \eqn{(y_{(i)},y_{(i+1)})}, the center is
#' \eqn{M=y_{(i)}+c(y_{(i+1)}-y_{(i)})} for the centrality parameter \code{c}.
#' For a given \eqn{c \in (0,1)}, the
#' expansion parameter \eqn{r} is taken to be \eqn{1/\max(c,1-c)} which yields non-degenerate asymptotic distribution of the
#' domination number.
#'
#' The test statistic is based on the binomial distribution, when success is defined as domination number being less than or
#' equal to 1 in the one interval case
#' (i.e., number of successes is equal to domination number \eqn{\le 1} in the partition intervals).
#' That is, the test statistic is based on the domination number for \code{Xp} points inside range of \code{Yp} points
#' (the domination numbers are summed over the \eqn{|Yp|-1} middle intervals)
#' for the PE-PCD and default end interval correction, \code{end.int.cor}, is \code{FALSE}
#' and the center \eqn{Mc} is chosen so that asymptotic distribution for the domination number is nondegenerate.
#' For this test to work, \code{Xp} must be at least 5 times more than \code{Yp} points
#' (or \code{Xp} must be at least 5 or more per partition interval).
#' Probability of success is the exact probability of success for the binomial distribution.
#'
#' **Caveat:** This test is currently a conditional test, where \code{Xp} points are assumed to be random, while \code{Yp} points are
#' assumed to be fixed (i.e., the test is conditional on \code{Yp} points).
#' Furthermore, the test is a large sample test when \code{Xp} points are substantially larger than \code{Yp} points,
#' say at least 7 times more.
#' This test is more appropriate when supports of \code{Xp} and \code{Yp} have a substantial overlap.
#' Currently, the \code{Xp} points outside the range of \code{Yp} points are handled with an end interval correction factor
#' (see the description below and the function code.)
#' Removing the conditioning and extending it to the case of non-concurring supports is
#' an ongoing line of research of the author of the package.
#'
#' See also (\insertCite{ceyhan:stat-2020;textual}{pcds}) for more on the uniformity test based on the arc
#' density of PE-PCDs.
#'
#' @param Xp A set of 1D points which constitute the vertices of the PE-PCD.
#' @param Yp A set of 1D points which constitute the end points of the partition intervals.
#' @param support.int Support interval \eqn{(a,b)} with \eqn{a<b}. Uniformity of \code{Xp} points in this interval
#' is tested.
#' @param c A positive real number which serves as the centrality parameter in PE proximity region;
#' must be in \eqn{(0,1)} (default \code{c=.5}).
#' @param end.int.cor A logical argument for end interval correction, default is \code{FALSE},
#' recommended when both \code{Xp} and \code{Yp} have the same interval support.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{\le 1)} for PE-PCD whose vertices are the 1D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding \code{alternative}.}
#' \item{conf.int}{Confidence interval for \eqn{Pr(}domination number\eqn{\le 1)} at the given level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{A \code{vector} with two entries: first is is the estimate of the parameter, i.e.,
#' \eqn{Pr(}domination number\eqn{\le 1)} and second is the domination number}
#' \item{null.value}{Hypothesized value for the parameter, i.e., the null value for \eqn{Pr(}domination number\eqn{\le 1)}}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{TSDomPEBin}} and \code{\link{PEdom1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10; supp<-c(a,b)
#' c<-.4
#'
#' r<-1/max(c,1-c)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#' TSDomPEBin1D(Xp,Yp,supp,c)
#' TSDomPEBin1D(Xp,Yp,supp,c,alt="l")
#' TSDomPEBin1D(Xp,Yp,supp,c,alt="g")
#' TSDomPEBin1D(Xp,Yp,supp,c,end=TRUE)
#' }
#'
#' @export TSDomPEBin1D
TSDomPEBin1D <- function(Xp,Yp,support.int,c=.5,end.int.cor=FALSE,
                           alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

  if ((!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp))))
  {stop('Xp and Yp must be 1D vectors of numerical entries.')}

  if (length(Yp)<2)
  {stop('Yp must be of length > 2')}

  if (!is.point(support.int) || support.int[2]<=support.int[1])
  {stop('support.int must be an interval as (a,b) with a<b')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  rstar<-1/max(c,1-c)  #r value for the non-degenerate asymptotic distribution

  p<-1-PG2PE1D.asy(c)  #asymptotic probability of success

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  nx<-length(Xp)  #number of Xp points
  ny<-length(Yp)  #number of Yp points
  nint<-ny-1  #number of middle intervals
  Ys<-sort(Yp)  #sorted Yp points (ends of the partition intervals)

  dom.num = PEdom1D(Xp,Yp,rstar,c)
  Gammas<-dom.num$int #domination numbers for the partition intervals
  #nint=ny-1 #-sum(Gammas0<1)

  Gam.all<-dom.num$d #domination number (with the end intervals included)
  Gam<-Gam.all-sum(sum(Xp<Ys[1])>0)-sum(sum(Xp>Ys[ny])>0)  #removing the domination number in the end intervals
  estimate2<-Gam

  estimate1<-Gam.all #domination number of the entire PE-PCD
  # ind0<- Gammas0>0 ; Gammas=Gammas0[ind0]
  #Bm.all<-sum(Gammas<=1) #the binomial test statistic, probability of success is Gamma <=1 for all intervals
  Bm<-sum(Gammas[-c(1,ny+1)]<=1) #the binomial test statistic, probability of success is Gamma <=1 for middle intervals

  method <-c("Large Sample Binomial Test based on the Domination Number of PE-PCD for Testing Uniformity of 1D Data ---")

  if (end.int.cor==TRUE)  #the part for the end interval correction
  {
    out.int<-sum(Xp<Ys[1])+sum(Xp>Ys[ny])
    prop.out<-out.int/nx #observed proportion of points in the end intervals
    exp.prop.out<-2/(ny+1)  #expected proportion of points in the end intervals

    x<-round(Bm*(1-(prop.out-exp.prop.out)))
    method <-c(method, " with End Interval Correction")
  } else
  { method <-c(method, " without End Interval Correction")}

  x<-Bm
  pval <-switch(alternative, less = pbinom(x, nint, p),
                greater = pbinom(x - 1, nint, p, lower.tail = FALSE),
                two.sided = {if (p == 0) (x == 0) else if (p == 1) (x == nint)
                  else { relErr <-1 + 1e-07
                  d <-dbinom(x, nint, p)
                  m <-nint * p
                  if (x == m) 1 else if (x < m)
                  {i <-seq.int(from = ceiling(m), to = nint)
                  y <-sum(dbinom(i, nint, p) <= d * relErr)
                  pbinom(x, nint, p) + pbinom(nint - y, nint, p, lower.tail = FALSE)
                  } else {
                    i <-seq.int(from = 0, to = floor(m))
                    y <-sum(dbinom(i, nint, p) <= d * relErr)
                    pbinom(y - 1, nint, p) + pbinom(x - 1, nint, p, lower.tail = FALSE)
                  }
                  }
                })

  p.L <- function(x, alpha) {
    if (x == 0)
      0
    else qbeta(alpha, x, nint - x + 1)
  }
  p.U <- function(x, alpha) {
    if (x == nint)
      1
    else qbeta(1 - alpha, x + 1, nint - x)
  }
  cint <-switch(alternative, less = c(0, p.U(x, 1 - conf.level)),
                greater = c(p.L(x, 1 - conf.level), 1), two.sided = {
                  alpha <-(1 - conf.level)/2
                  c(p.L(x, alpha), p.U(x, alpha))
                })
  attr(cint, "conf.level") <- conf.level

  estimate2 <-x/nint

  #names(x) <-ifelse(end.int.cor==TRUE,"corrected domination number","domination number" )
  names(x) <- "#(domination number is <= 1)" #"domination number - number of partition intervals"
  names(nint) <-"number of partition (middle) intervals based on Yp"
  names(p) <-"Pr(Domination Number <= 1)"

  names(estimate1) <-c(" domination number ")
  names(estimate2) <-c("|| Pr(domination number <= 1)")

  structure(
    list(statistic = x,
         p.value = pval,
         conf.int = cint,
         estimate = c(estimate1,estimate2),
         null.value = p,
         alternative = alternative,
         method = method,
         data.name = dname
    ),
    class = "htest")

} #end of the function
