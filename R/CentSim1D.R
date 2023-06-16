##################################################################################################################
################################CENTRAL SIMILARITY FUNCTIONS######################################################
##################################################################################################################
#CentSim1D.R;
#################################################################
#Functions for NCS in R^1
#################################################################

#' @title The indicator for the presence of an arc from a point to another for Central Similarity
#' Proximity Catch Digraphs (CS-PCDs) - middle interval case
#'
#' @description Returns \eqn{I(p_2} in \eqn{N_{CS}(p_1,t,c))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{CS}(p_1,t,c)}, returns 0
#' otherwise, where \eqn{N_{CS}(x,t,c)} is the CS proximity region for point \eqn{x} and is constructed with expansion
#' parameter \eqn{t>0} and centrality parameter \eqn{c \in (0,1)} for the interval \eqn{(a,b)}.
#'
#' CS proximity regions are defined with respect to the middle interval \code{int} and vertex regions are based
#' on the center associated with the centrality parameter \eqn{c \in (0,1)}. For the interval, \code{int}\eqn{=(a,b)}, the
#' parameterized center is \eqn{M_c=a+c(b-a)}. \code{rv} is the index of the vertex region \eqn{p_1} resides, with default=\code{NULL}.
#'
#' If \eqn{p_1} and \eqn{p_2} are distinct and either of them are outside interval \code{int}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations
#' (i.e., loops are allowed in the digraph).
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param p1,p2 1D points; \eqn{p_1} is the point for which the proximity region, \eqn{N_{CS}(p_1,t,c)} is
#' constructed and \eqn{p_2} is the point which the function is checking whether its inside
#' \eqn{N_{CS}(p_1,t,c)} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param rv Index of the end interval containing the point, either \code{1,2} or \code{NULL} (default is \code{NULL}).
#'
#' @return \eqn{I(p_2} in \eqn{N_{CS}(p_1,t,c))} for points \eqn{p_1} and \eqn{p_2} that is, returns 1 if \eqn{p_2} is in \eqn{N_{CS}(p_1,t,c)},
#' returns 0 otherwise
#'
#' @seealso \code{\link{IarcCSend.int}}, \code{\link{IarcPEmid.int}}, and \code{\link{IarcPEend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.5
#' t<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' IarcCSmid.int(7,5,int,t,c)
#' IarcCSmid.int(7,7,int,t,c)
#' IarcCSmid.int(7,5,int,t,c=.4)
#'
#' IarcCSmid.int(1,3,int,t,c)
#'
#' IarcCSmid.int(9,11,int,t,c)
#'
#' IarcCSmid.int(19,1,int,t,c)
#' IarcCSmid.int(19,19,int,t,c)
#'
#' IarcCSmid.int(3,5,int,t,c)
#'
#' #or try
#' Rv<-rel.vert.mid.int(3,int,c)$rv
#' IarcCSmid.int(3,5,int,t,c,rv=Rv)
#'
#' IarcCSmid.int(7,5,int,t,c)
#'
#' @export IarcCSmid.int
IarcCSmid.int <- function(p1,p2,int,t,c=.5,rv=NULL)
{
  if (!is.point(p1,1) || !is.point(p2,1) )
  {stop('p1 and p2 must be scalars')}

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (p1==p2 )
  {arc<-1; return(arc); stop}

  y1<-int[1]; y2<-int[2];
  if (p1<y1 || p1>y2 || p2<y1 || p2>y2 )
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  {rv<-rel.vert.mid.int(p1,int,c)$rv #determines the vertex region for 1D point p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2))!=1)
  {stop('vertex index, rv, must be 1 or 2')}}

  arc<-0;
  if (rv==1)
  {
    if ( p2 < p1+t*(1-c)*(p1-y1)/c & p2 > p1-t*(p1-y1) ) {arc <-1}
  } else {
    if ( p2 < p1+t*(y2-p1) & p2 > p1-c*t*(y2-p1)/(1-c) ) {arc<-1}
  }
  arc
} #end of the function
#'

#################################################################

#' @title Number of Arcs of of Central Similarity Proximity Catch Digraphs (CS-PCDs) - middle interval case
#'
#' @description Returns the number of arcs of of Central Similarity Proximity Catch Digraphs (CS-PCDs) whose vertices are
#' the given 1D numerical data set, \code{Xp}.
#'
#' CS proximity region \eqn{N_{CS}(x,t,c)} is defined with respect to the interval
#' \code{int}\eqn{=(a,b)} for this function. CS proximity region is constructed with expansion parameter \eqn{t>0} and
#' centrality parameter \eqn{c \in (0,1)}.
#'
#' Vertex regions are based on the center associated with the centrality parameter \eqn{c \in (0,1)}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)} and for the number of arcs,
#' loops are not allowed so arcs are only possible for points inside the middle interval \code{int} for this function.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return Number of arcs for the PE-PCD whose vertices are the 1D data set, \code{Xp},
#' with expansion parameter, \eqn{r \ge 1}, and centrality parameter, \eqn{c \in (0,1)}. PE proximity regions are defined only
#' for \code{Xp} points inside the interval \code{int}, i.e., arcs are possible for such points only.
#'
#' @seealso \code{\link{num.arcsCSend.int}}, \code{\link{num.arcsPEmid.int}}, and \code{\link{num.arcsPEend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' t<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-10
#' Xp<-runif(n,a,b)
#' num.arcsCSmid.int(Xp,int,t,c)
#'
#' num.arcsCSmid.int(Xp,int,t,c=.3)
#'
#' num.arcsCSmid.int(Xp,int,t=1.5,c)
#'
#' #num.arcsCSmid.int(Xp,int,t,c+5) #gives error
#' #num.arcsCSmid.int(Xp,int,t,c+10)
#'
#' n<-10  #try also n<-20
#' Xp<-runif(n,a-5,b+5)
#' num.arcsCSint(Xp,int,t,c)
#'
#' Xp<-runif(n,a+10,b+10)
#' num.arcsCSmid.int(Xp,int,t,c)
#'
#' n<-10
#' Xp<-runif(n,a,b)
#' num.arcsCSmid.int(Xp,int,t,c)
#'
#' @export num.arcsCSmid.int
num.arcsCSmid.int <- function(Xp,int,t,c=.5)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp vector of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  n<-length(Xp)
  if (n>0)
  {
    arcs<-0
    for (i in 1:n)
    {x1<-Xp[i]
    if (x1>y2 || x1<y1)
    {
      arcs<-arcs+0
    } else
    {
      v<-rel.vert.mid.int(x1,int,c)$rv
      if (v==1)
      {
        xR<-x1+t*(1-c)*(x1-y1)/c
        xL<-x1-t*(x1-y1)
        arcs<-arcs+sum((Xp < min(xR,y2)) & (Xp > max(xL,y1)))-1 #minus 1 is for loops
      } else {
        xR <-x1+t*(y2-x1)
        xL <-x1-c*t*(y2-x1)/(1-c)
        arcs<-arcs+sum((Xp < min(xR,y2)) & (Xp > max(xL,y1)))-1
      }
    }
    }
  } else
  {arcs<-0}
  arcs
} #end of the function
#'

#################################################################

# funsMuVarCS1D
#'
#' @title Returning the mean and (asymptotic) variance of arc density of Central Similarity Proximity Catch Digraph (CS-PCD)
#' for 1D data - middle interval case
#'
#' @description
#' Two functions: \code{muCS1D} and \code{asyvarCS1D}.
#'
#' \code{muCS1D} returns the mean of the (arc) density of CS-PCD
#' and \code{asyvarCS1D} returns the (asymptotic) variance of the arc density of CS-PCD
#' for a given centrality parameter \eqn{c \in (0,1)} and an expansion parameter \eqn{t>0} and 1D uniform data in a
#' finite interval \eqn{(a,b)}, i.e., data from \eqn{U(a,b)} distribution.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return \code{muCS1D} returns the mean and \code{asyvarCS1D} returns the asymptotic variance of the
#' arc density of CS-PCD for uniform data in an interval
#'
#' @name funsMuVarCS1D
NULL
#'
#' @seealso \code{\link{muPE1D}} and \code{\link{asyvarPE1D}}
#'
#' @rdname funsMuVarCS1D
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \donttest{
#' #Examples for muCS1D
#' muCS1D(1.2,.4)
#' muCS1D(1.2,.6)
#'
#' tseq<-seq(0.01,5,by=.05)
#' cseq<-seq(0.01,.99,by=.05)
#'
#' ltseq<-length(tseq)
#' lcseq<-length(cseq)
#'
#' mu.grid<-matrix(0,nrow=ltseq,ncol=lcseq)
#' for (i in 1:ltseq)
#'   for (j in 1:lcseq)
#'   {
#'     mu.grid[i,j]<-muCS1D(tseq[i],cseq[j])
#'   }
#'
#' persp(tseq,cseq,mu.grid, xlab="t", ylab="c", zlab="mu(t,c)",theta = -30,
#' phi = 30, expand = 0.5, col = "lightblue", ltheta = 120,
#' shade = 0.05, ticktype = "detailed")
#' }
#'
#' @export
muCS1D <- function(t,c)
{
  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  mean<-0;
  if (t<1)
  {
    mean<-t/2;
  } else {
    mean<-1/2*t*(2*c^2*t-2*c^2-2*c*t+2*c-1)/((c*t-c-t)*(c*t-c+1));
  }
  mean
} #end of the function
#'
#' @rdname funsMuVarCS1D
#'
#' @examples
#' \donttest{
#' #Examples for asyvarCS1D
#' asyvarCS1D(1.2,.8)
#'
#' tseq<-seq(0.01,5,by=.05)
#' cseq<-seq(0.01,.99,by=.05)
#'
#' ltseq<-length(tseq)
#' lcseq<-length(cseq)
#'
#' var.grid<-matrix(0,nrow=ltseq,ncol=lcseq)
#' for (i in 1:ltseq)
#'   for (j in 1:lcseq)
#'   {
#'     var.grid[i,j]<-asyvarCS1D(tseq[i],cseq[j])
#'   }
#'
#' persp(tseq,cseq,var.grid, xlab="t", ylab="c", zlab="var(t,c)", theta = -30,
#' phi = 30, expand = 0.5, col = "lightblue", ltheta = 120,
#' shade = 0.05, ticktype = "detailed")
#' }
#'
#' @export
asyvarCS1D <- function(t,c)
{
  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  asyvar<-0;
  if (t<1)
  {
    asyvar<--1/3*t^2*(c^2*t^3-3*c^2*t^2-c*t^3+2*c^2*t+3*c*t^2-c^2-2*c*t-t^2+c+t)/((c*t-c+1)*(c*t-c-t));
  } else {
    asyvar<--1/3*c*(c-1)*(2*c^4*t^5-7*c^4*t^4-4*c^3*t^5+8*c^4*t^3+14*c^3*t^4+3*c^2*t^5-2*c^4*t^2-16*c^3*t^3-7*c^2*t^4-c*t^5-2*c^4*t+4*c^3*t^2+12*c^2*t^3+c^4+4*c^3*t-6*c^2*t^2-4*c*t^3-2*c^3-3*c^2*t+4*c*t^2+c^2+c*t-t^2)/((c*t-c-t)^3*(c*t-c+1)^3);
  }
  asyvar
} #end of the function
#'

#################################################################

#' @title A test of uniformity of 1D data in a given interval based on Central Similarity Proximity Catch Digraph
#' (CS-PCD)
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of uniformity of 1D
#' data in one interval based on the normal approximation of the arc density of the CS-PCD with expansion
#' parameter \eqn{t>0} and centrality parameter \eqn{c \in (0,1)}.
#'
#' The function yields the test statistic, \eqn{p}-value for the
#' corresponding \code{alternative}, the confidence interval, estimate and null value for the parameter of interest
#' (which is the arc density), and method and name of the data set used.
#'
#' The null hypothesis is that data is
#' uniform in a finite interval (i.e., arc density of CS-PCD equals to its expected value under uniform
#' distribution) and \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the end
#' points) or right-sided (i.e., data is accumulated around the mid point or center \eqn{M_c}).
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the arc density of CS-PCD based on
#' the 1D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval for the arc density at the given level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate of the parameter, i.e., arc density}
#' \item{null.value}{Hypothesized value for the parameter, i.e., the null arc density, which is usually the
#' mean arc density under uniform distribution.}
#' \item{alternative}{Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{PEarc.dens.test.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \donttest{
#' c<-.4
#' t<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-10
#' Xp<-runif(n,a,b)
#'
#' num.arcsCSmid.int(Xp,int,t,c)
#' CSarc.dens.test.int(Xp,int,t,c)
#'
#' num.arcsCSmid.int(Xp,int,t,c=.3)
#' CSarc.dens.test.int(Xp,int,t,c=.3)
#'
#' num.arcsCSmid.int(Xp,int,t=1.5,c)
#' CSarc.dens.test.int(Xp,int,t=1.5,c)
#'
#' Xp<-runif(n,a-1,b+1)
#' num.arcsCSmid.int(Xp,int,t,c)
#' CSarc.dens.test.int(Xp,int,t,c)
#'
#' c<-.4
#' t<-.5
#' a<-0; b<-10; int<-c(a,b)
#' n<-10  #try also n<-20
#' Xp<-runif(n,a,b)
#'
#' CSarc.dens.test.int(Xp,int,t,c)
#' }
#'
#' @export CSarc.dens.test.int
CSarc.dens.test.int <- function(Xp,int,t,c=.5,alternative = c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

  if (!is.point(Xp,length(Xp)))
  {stop('Xp vector of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
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

  num.arcs<-num.arcsCSmid.int(Xp,int,t,c)
  arc.dens<-num.arcs/(n*(n-1))
  estimate1<-arc.dens
  mn<-muCS1D(t,c)
  asy.var<-asyvarCS1D(t,c)
  TS<-sqrt(n) *(arc.dens-mn)/sqrt(asy.var)
  method <-c("Large Sample z-Test Based on Arc Density of CS-PCD for Testing Uniformity of 1D Data")
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

  attr(cint, "conf.level") <- conf.level

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

#' @title A test of segregation/association based on arc density of Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 1D data
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points in the range (i.e., range) of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster away from \code{Yp} points) and association (where \code{Xp} points cluster around
#' \code{Yp} points) based on the normal approximation of the arc density of the CS-PCD for uniform 1D data.
#'
#' The function yields the test statistic, \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval, estimate and null value for the parameter of interest (which is the arc density),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points in the range of \code{Yp} points, arc density
#' of CS-PCD whose vertices are \code{Xp} points equals to its expected value under the uniform distribution and
#' \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated around the centers of the intervals, or segregation).
#'
#' CS proximity region is constructed with the expansion parameter \eqn{t > 0} and centrality parameter \code{c} which yields
#' \eqn{M}-vertex regions. More precisely, for a middle interval \eqn{(y_{(i)},y_{(i+1)})}, the center is
#' \eqn{M=y_{(i)}+c(y_{(i+1)}-y_{(i)})} for the centrality parameter \eqn{c \in (0,1)}.
#' This test is more appropriate when supports of \code{Xp} and \code{Yp} has a substantial overlap.
#'
#' \code{end.int.cor} is for end interval correction, (default is "no end interval correction", i.e., \code{end.int.cor=FALSE}),
#' recommended when both \code{Xp} and \code{Yp} have the same interval support.
#'
#' @param Xp A set of 1D points which constitute the vertices of the CS-PCD.
#' @param Yp A set of 1D points which constitute the end points of the partition intervals.
#' @param support.int Support interval \eqn{(a,b)} with \eqn{a<b}.
#' Uniformity of \code{Xp} points in this interval is tested. Default is \code{NULL}.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number which serves as the centrality parameter in CS proximity region;
#' must be in \eqn{(0,1)} (default \code{c=.5}).
#' @param end.int.cor A logical argument for end interval correction, default is \code{FALSE},
#' recommended when both \code{Xp} and \code{Yp} have the same interval support.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the arc density
#' CS-PCD whose vertices are the 1D data set \code{Xp}.
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
#' @seealso \code{\link{CSarc.dens.test}} and \code{\link{CSarc.dens.test.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' tau<-2
#' c<-.4
#' a<-0; b<-10; int=c(a,b)
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
#' CSarc.dens.test1D(Xp,Yp,tau,c,int)
#' CSarc.dens.test1D(Xp,Yp,tau,c,int,alt="l")
#' CSarc.dens.test1D(Xp,Yp,tau,c,int,alt="g")
#'
#' CSarc.dens.test1D(Xp,Yp,tau,c,int,end.int.cor = TRUE)
#'
#' Yp2<-runif(ny,a,b)+11
#' CSarc.dens.test1D(Xp,Yp2,tau,c,int)
#'
#' n<-10  #try also n<-20
#' Xp<-runif(n,a,b)
#' CSarc.dens.test1D(Xp,Yp,tau,c,int)
#'
#' @export CSarc.dens.test1D
CSarc.dens.test1D <- function(Xp,Yp,t,c=.5,support.int=NULL,end.int.cor=FALSE,
                            alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

  if ((!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp))))
  {stop('Xp and Yp must be 1D vectors of numerical entries.')}

  nx =length(Xp); ny = length(Yp)
  if (ny<2)
  {stop('Yp must be of length > 2')}

  if (!is.null(support.int))
  {
  if (!is.point(support.int) || support.int[2]<=support.int[1])
  {stop('support.int must be an interval as (a,b) with a<b')}
  }

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(t,1) || t <=0)
  {stop('t must be a scalar greater than 0')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  Arcs<-num.arcsCS1D(Xp,Yp,t,c)  #uses the default c=.5, unless specified otherwise
  NinR<-Arcs$num.in.range #number of Xp points in the range of Yp points

  num.arcs.ints = Arcs$int.num.arcs #vector of number of arcs in the partition intervals
  n.int = length(num.arcs.ints)
  num.arcs = sum(num.arcs.ints[-c(1,n.int)]) #this is to remove the number of arcs in the end intervals
  num.dat.ints = Arcs$num.in.int[-c(1,n.int)] #vector of numbers of data points in the partition intervals
  Wvec<-Arcs$w
  LW<-Wvec/sum(Wvec)

  dat.int.ind = Arcs$data.int.ind #indices of partition intervals in which data points reside
  mid.ind = which(dat.int.ind!=1 & dat.int.ind!=n.int) #indices of Xp points in range of Yp points (i.e., in middle intervals)
  dat.int.ind = dat.int.ind[mid.ind] #removing the end interval indices
  dat.mid = Xp[mid.ind] #Xp points in range of Yp points (i.e., in middle intervals)
  part.int.mid =  t(Arcs$partition.int)[-c(1,n.int),] #middle partition intervals

  ind.Xp1 = which(num.dat.ints==1) #indices of partition intervals containing only one Xp point
  if (length(ind.Xp1)>0)
  {
    for (i in ind.Xp1)
    {
      Xpi = dat.mid[dat.int.ind==i+1]
      int =  part.int.mid[i,]
      npe = NCSint(Xpi,int,t,c)
      num.arcs = num.arcs+(npe[2]-npe[1])/Wvec[i]
    }
  }

  asy.mean0<-muCS1D(t,c)  #asy mean value for the (t,c) pair
  asy.mean<-asy.mean0*sum(LW^2)

  asy.var0<-asyvarCS1D(t,c)  #asy variance value for the (t,c) pair
  asy.var<-asy.var0*sum(LW^3)+4*asy.mean0^2*(sum(LW^3)-(sum(LW^2))^2)

  if (NinR  == 0)
  {warning('There is no Xp point in the range of Yp points to compute arc density,
           but as this is clearly a segregation pattern, arc density is taken to be 1!')
    arc.dens=1
    TS0<-sqrt(nx)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat
  } else
  {  arc.dens<-num.arcs/(NinR*(NinR-1))
  TS0<-sqrt(NinR)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat}  #arc density
  }
  estimate1<-arc.dens
  estimate2<-asy.mean

  method <-c("Large Sample z-Test Based on Arc Density of CS-PCD for Testing Uniformity of 1D Data ---")
  if (end.int.cor==F)
  {
    TS<-TS0
    method <-c(method,"\n without End Interval Correction")
  }
  else
  {
    n<-length(Xp)  #number of X points
    m<-length(Yp)  #number of Y points
    NoutRange<-n-NinR #number of points outside of the range

    prop.out<-NoutRange/n #observed proportion of points outside range
    exp.prop.out<-2/m  #expected proportion of points outside range of Y points

    TS<-TS0+abs(TS0)*sign(prop.out-exp.prop.out)*(prop.out-exp.prop.out)^2
    method <-c(method,"\n with End Interval Correction")
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

  attr(cint, "conf.level") <- conf.level

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

#' @title The arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for 1D data - middle intervals case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of CS-PCD and related parameters and the quantities of the digraph.
#'
#' For this function, CS proximity regions are constructed with respect to the intervals
#' based on \code{Yp} points with expansion parameter \eqn{t>0} and centrality parameter \eqn{c \in (0,1)}. That is, for this
#' function, arcs may exist for points only inside the intervals.
#' It also provides various descriptions and quantities about the arcs of the CS-PCD
#' such as number of arcs, arc density, etc.
#'
#' Vertex regions are based on center \eqn{M_c} of each middle interval.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the CS-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the intervals.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, they are expansion and centrality parameters.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, i.e., \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of CS-PCD for 1D data in the middle intervals}
#' \item{E}{Heads (or arrow ends) of the arcs of CS-PCD for 1D data in the middle intervals}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsPEend.int}}, \code{\link{arcsPE1D}}, \code{\link{arcsCSmid.int}},
#'  \code{\link{arcsCSend.int}} and \code{\link{arcsCS1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-1.5
#' c<-.4
#' a<-0; b<-10
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' arcsCSmid.int(Xp,Yp,t,c)
#' arcsCSmid.int(Xp,Yp+10,t,c)
#'
#' Arcs<-arcsCSmid.int(Xp,Yp,t,c)
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
#' Xlim<-range(Xp,Yp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),
#' main="arcs of CS-PCD whose vertices (jittered along y-axis)\n in middle intervals ",
#' xlab=" ", ylab=" ", xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit),pch=".")
#' abline(h=0,lty=1)
#' points(Xp, yjit,pch=".",cex=3)
#' abline(v=Yp,lty=2)
#' arrows(S, yjit, E, yjit, length = .05, col= 4)
#'
#' t<-.5
#' c<-.4
#' a<-0; b<-10;
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' arcsCSmid.int(Xp,Yp,t,c)
#'
#' @export arcsCSmid.int
arcsCSmid.int <- function(Xp,Yp,t,c=.5)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
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

      #the arcs of CS-PCDs for parameters t and c
      S<-E<-vector()  #S is for source and E is for end points for the arcs for middle intervals
      for (i in 1:nt)
      {
        Xi<-Xint[i.int==i] #X points in the ith Yp mid interval
        ni<-length(Xi)
        if (ni>1 )
        {
          y1<-Ys[i]; y2<-Ys[i+1]; int<-c(y1,y2)
          for (j in 1:ni)
          {x1 <-Xi[j]; Xinl<-Xi[-j] #to avoid loops
          v<-rel.vert.mid.int(x1,int,c)$rv
          if (v==1)
          {
            xR<-x1+t*(1-c)*(x1-y1)/c
            xL<-x1-t*(x1-y1)
            ind.tails<-((Xinl < min(xR,y2)) & (Xinl > max(xL,y1)))
            st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
            S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
          } else {
            xR <-x1+t*(y2-x1)
            xL <-x1-c*t*(y2-x1)/(1-c)
            ind.tails<-((Xinl < min(xR,y2)) & (Xinl > max(xL,y1)))
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

  param<-c(t,c)
  names(param)<-c("expansion parameter","centrality parameter")
  typ<-paste("Central Similarity Proximity Catch Digraph (CS-PCD) for 1D Points in the Middle Intervals with Expansion Parameter t = ",t," and Centrality Parameter c = ",c,sep="")

  main.txt<-paste("Arcs of CS-PCD for Points (jittered\n along y-axis) in Middle Intervals with t = ",round(t,2)," and c = ",round(c,2),sep="")

  nvert<-nx2; nint<-ny-1; narcs=ifelse(!is.null(S),length(S),0);
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
#The case of end intervals
#################################################################

#' @title The indicator for the presence of an arc from a point to another for
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) - end interval case
#'
#' @description Returns \eqn{I(p_2} in \eqn{N_{CS}(p_1,t))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{CS}(p_1,t)}, returns 0
#' otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{t>0}
#' for the region outside the interval \eqn{(a,b)}.
#'
#' \code{rv} is the index of the end vertex region \eqn{p_1} resides, with default=\code{NULL},
#' and \code{rv=1} for left end interval and \code{rv=2} for the right end interval.
#' If \eqn{p_1} and \eqn{p_2} are distinct and either of them are inside interval \code{int}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param p1 A 1D point for which the CS proximity region is constructed.
#' @param p2 A 1D point to check whether it is inside the proximity region or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param rv Index of the end interval containing the point, either \code{1,2} or \code{NULL} (default=\code{NULL}).
#'
#' @return \eqn{I(p_2} in \eqn{N_{CS}(p_1,t))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{CS}(p_1,t)}
#' (i.e., if there is an arc from \eqn{p_1} to \eqn{p_2}), returns 0 otherwise
#'
#' @seealso \code{\link{IarcCSmid.int}}, \code{\link{IarcPEmid.int}}, and \code{\link{IarcPEend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' a<-0; b<-10; int<-c(a,b)
#' t<-2
#'
#' IarcCSend.int(15,17,int,t)
#' IarcCSend.int(15,15,int,t)
#'
#' IarcCSend.int(1.5,17,int,t)
#' IarcCSend.int(1.5,1.5,int,t)
#'
#' IarcCSend.int(-15,17,int,t)
#'
#' IarcCSend.int(-15,-17,int,t)
#'
#' a<-0; b<-10; int<-c(a,b)
#' t<-.5
#'
#' IarcCSend.int(15,17,int,t)
#'
#' @export IarcCSend.int
IarcCSend.int <- function(p1,p2,int,t,rv=NULL)
{
  if (!is.point(p1,1) || !is.point(p2,1) )
  {stop('p1 and p2 must be scalars')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (p1==p2 )
  {arc<-1; return(arc); stop}

  if ((p1>y1 & p1<y2) || (p2>y1 & p2<y2))
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  {rv<-rel.vert.end.int(p1,int)$rv #determines the vertex for the end interval for 1D point p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2))!=1)
  {stop('vertex index, rv, must be 1 or 2')}}

  arc<-0;
  if (rv==1)
  {
    if ( p2 < p1+t*(y1-p1) & p2 > p1-t*(y1-p1) ) {arc <-1}
  } else {
    if ( p2 < p1+t*(p1-y2) & p2 > p1-t*(p1-y2) ) {arc<-1}
  }
  arc
} #end of the function
#'

#################################################################

#' @title Number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs) - end interval case
#'
#' @description Returns the number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs) whose
#' vertices are a 1D numerical data set, \code{Xp}, outside the interval \code{int}\eqn{=(a,b)}.
#'
#' CS proximity region is constructed only with expansion parameter \eqn{t>0} for points outside the interval \eqn{(a,b)}.
#'
#' End vertex regions are based on the end points of the interval,
#' i.e., the corresponding end vertex region is an interval as \eqn{(-\infty,a)} or \eqn{(b,\infty)} for the interval \eqn{(a,b)}.
#' For the number of arcs, loops are not allowed, so arcs are only possible for points outside
#' the interval, \code{int}, for this function.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A \code{vector} of 1D points which constitute the vertices of the digraph.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return Number of arcs for the CS-PCD with vertices being 1D data set, \code{Xp},
#' expansion parameter, \code{t}, for the end intervals.
#'
#' @seealso \code{\link{num.arcsCSmid.int}}, \code{\link{num.arcsPEmid.int}}, and \code{\link{num.arcsPEend.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-5
#' XpL<-runif(n,a-5,a)
#' XpR<-runif(n,b,b+5)
#' Xp<-c(XpL,XpR)
#'
#' num.arcsCSend.int(Xp,int,t=2)
#'
#' num.arcsCSend.int(Xp,int,t=1.2)
#'
#' num.arcsCSend.int(Xp,int,t=4)
#'
#' num.arcsCSend.int(Xp,int,t=2+5)
#' #num.arcsCSend.int(Xp,int,t=c(-5,15))
#'
#' n<-10  #try also n<-20
#' Xp2<-runif(n,a-5,b+5)
#' num.arcsCSend.int(Xp2,int,t=2)
#'
#' t<-.5
#' num.arcsCSend.int(Xp,int,t)
#'
#' @export num.arcsCSend.int
num.arcsCSend.int <- function(Xp,int,t)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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
    {p1<-Xp[i]; rv<-rel.vert.end.int(p1,int)$rv
    for (j in ((1:n)[-i]) )
    {p2<-Xp[j]
    arcs<-arcs+IarcCSend.int(p1,p2,int,t,rv)
    }
    }
  }
  arcs
} #end of the function
#'

#################################################################

#' funsMuVarCSend.int
#'
#' @title Returns the mean and (asymptotic) variance of arc density of Central Similarity Proximity Catch Digraph (CS-PCD)
#' for 1D data - end interval case
#'
#' @description
#' Two functions: \code{muCSend.int} and \code{asyvarCSend.int}.
#'
#' \code{muCSend.int} returns the mean of the arc density of CS-PCD
#' and \code{asyvarCSend.int} returns the asymptotic variance of the arc density of CS-PCD
#' for a given expansion parameter \eqn{t>0}  for 1D uniform data in the left and right end intervals
#' for the interval \eqn{(a,b)}.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#'
#' @return \code{muCSend.int} returns the mean and \code{asyvarCSend.int} returns the asymptotic variance of the
#' arc density of CS-PCD for uniform data in end intervals
#'
#' @name funsMuVarCSend.int
NULL
#'
#' @seealso \code{\link{muPEend.int}} and \code{\link{asyvarPEend.int}}
#'
#' @rdname funsMuVarCSend.int
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for muCSend.int
#' muCSend.int(1.2)
#'
#' tseq<-seq(0.01,5,by=.05)
#' ltseq<-length(tseq)
#'
#' mu.end<-vector()
#' for (i in 1:ltseq)
#' {
#'   mu.end<-c(mu.end,muCSend.int(tseq[i]))
#' }
#'
#' oldpar <- par(no.readonly = TRUE)
#' par(mar = c(5,4,4,2) + 0.1)
#' plot(tseq, mu.end,type="l",
#' ylab=expression(paste(mu,"(t)")),xlab="t",lty=1,xlim=range(tseq),ylim=c(0,1))
#' par(oldpar)
#'
#' @export muCSend.int
muCSend.int <- function(t)
{
  if (!is.point(t,1) || t<=0)
  {stop('the argument must be a scalar greater than 0')}

  mean<-0;
  if (t<1)
  {
    mean<-t*(t+2)/(2*(1+t));
  } else {
    mean<-(1+2*t)/(2*(1+t));
  }
  mean
} #end of the function
#'
#' @rdname funsMuVarCSend.int
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #Examples for asyvarCSend.int
#' asyvarCSend.int(1.2)
#'
#' tseq<-seq(.01,5,by=.05)
#' ltseq<-length(tseq)
#'
#' var.end<-vector()
#' for (i in 1:ltseq)
#' {
#'   var.end<-c(var.end,asyvarCSend.int(tseq[i]))
#' }
#'
#' oldpar <- par(no.readonly = TRUE)
#' par(mar=c(5,5,4,2))
#' plot(tseq, var.end,type="l",xlab="t",ylab=expression(paste(sigma^2,"(t)")),lty=1,xlim=range(tseq))
#' par(oldpar)
#'
#' @export asyvarCSend.int
asyvarCSend.int <- function(t)
{
  if (!is.point(t,1) || t<=0)
  {stop('the argument must be a scalar greater than 0')}

  asyvar<-0;
  if (t<1)
  {
    asyvar<-t^2*(2*t^4+4*t^3+t^2-4*t-4)/(3*(1+t)^3);
  } else {
    asyvar<-t^2/(3*(1+t)^3);
  }
  asyvar
} #end of the function
#'

#################################################################

#' @title The arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for 1D data - end interval case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of CS-PCD and related parameters and the quantities of the digraph.
#' \code{Yp} determines the end points of the end intervals.
#'
#' For this function, CS proximity regions are constructed data points outside the intervals based on
#' \code{Yp} points with expansion parameter \eqn{t>0}. That is, for this function,
#' arcs may exist for points only inside end intervals.
#' It also provides various descriptions and quantities about the arcs of the CS-PCD
#' such as number of arcs, arc density, etc.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the CS-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the intervals.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, it is the expansion parameter.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization based on \code{Yp}.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitutes the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of CS-PCD for 1D data in the end intervals}
#' \item{E}{Heads (or arrow ends) of the arcs of CS-PCD for 1D data in the end intervals}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals (which is 2 for end intervals), number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsCSmid.int}}, \code{\link{arcsCS1D}} , \code{\link{arcsPEmid.int}},
#' \code{\link{arcsPEend.int}}  and \code{\link{arcsPE1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-1.5
#' a<-0; b<-10; int<-c(a,b)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xr<-range(a,b)
#' xf<-(xr[2]-xr[1])*.5
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' arcsCSend.int(Xp,Yp,t)
#'
#' Arcs<-arcsCSend.int(Xp,Yp,t)
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
#' main="arcs of CS-PCD with vertices (jittered along y-axis)\n in end intervals ",
#'      xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit))
#' abline(h=0,lty=1)
#' points(Xp, yjit,pch=".",cex=3)
#' abline(v=Yp,lty=2)
#' arrows(S, yjit, E, yjit, length = .05, col= 4)
#'
#' arcsCSend.int(Xp,Yp,t)
#'
#' @export arcsCSend.int
arcsCSend.int <- function(Xp,Yp,t)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))
  tname <-deparse(substitute(t))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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
    xR<-x1+t*(ymin-x1); xL<-x1-t*(ymin-x1)
    ind.tails<-((xLe < min(xR,ymin)) & (xLe > xL))
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
    xR<-x1+t*(x1-ymax); xL<-x1-t*(x1-ymax)
    ind.tails<-((xRe < xR) & xRe > max(ymax,xL))
    st<-sum(ind.tails)  #sum of tails of the arcs with head XRe[j]
    S<-c(S,rep(x1,st)); E<-c(E,xRe[ind.tails])
    }
  }

  if (length(S)==0)
  {S<-E<-NA}

  param<-t
  names(param)<-"expansion parameter"
  typ<-paste("Central Similarity Proximity Catch Digraph (CS-PCD) for 1D Points in the End Intervals with Expansion Parameter t = ",t,sep="")

  main.txt<-paste("Arcs of CS-PCD for Points (jittered\n along y-axis) in End Intervals with t = ",round(t,2),sep="")

  nvert<-nle+nre; ny<-length(Yp); nint<-2; narcs=ifelse(!is.null(S),length(S),0);
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

#' @title The plot of the arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs) for 1D data
#' (vertices jittered along \eqn{y}-coordinate) - multiple interval case
#'
#' @description Plots the arcs of CS-PCD whose vertices are the 1D points, \code{Xp}. CS proximity regions are constructed with
#' expansion parameter \eqn{t>0} and centrality parameter \eqn{c \in (0,1)} and the intervals are based on \code{Yp} points (i.e.
#' the intervalization is based on \code{Yp} points). That is, data set \code{Xp}
#' constitutes the vertices of the digraph and \code{Yp} determines the end points of the intervals.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default for \eqn{Jit=.1}) is added to
#' the \eqn{y}-direction where \code{Jit} equals to the range of \code{Xp} and \code{Yp} multiplied by \code{Jit} with default for \eqn{Jit=.1}).
#'
#' \code{centers} is a logical argument, if \code{TRUE}, plot includes the centers of the intervals
#' as vertical lines in the plot, else centers of the intervals are not plotted.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A \code{vector} of 1D points constituting the vertices of the CS-PCD.
#' @param Yp A \code{vector} of 1D points constituting the end points of the intervals.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis where \code{Jit} equals to the range of \code{Xp} and \code{Yp} multiplied by
#' \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles of the \eqn{x} and \eqn{y} axes in the plot (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param centers A logical argument, if \code{TRUE}, plot includes the centers of the intervals
#' as vertical lines in the plot, else centers of the intervals are not plotted.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of CS-PCD whose vertices are the 1D data set \code{Xp} in which vertices are jittered
#' along \eqn{y}-axis for better visualization.
#'
#' @seealso \code{\link{plotPEarcs1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-1.5
#' c<-.4
#' a<-0; b<-10; int<-c(a,b)
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xr<-range(a,b)
#' xf<-(xr[2]-xr[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' Xlim=range(Xp,Yp)
#' Ylim=c(-.2,.2)
#'
#' jit<-.1
#'
#' plotCSarcs1D(Xp,Yp,t,c,jit,xlab="",ylab="",xlim=Xlim,ylim=Ylim)
#'
#' set.seed(1)
#' plotCSarcs1D(Xp,Yp,t=1.5,c=.3,jit,main="t=1.5, c=.3",xlab="",ylab="",centers=TRUE)
#' set.seed(1)
#' plotCSarcs1D(Xp,Yp,t=2,c=.3,jit,main="t=2, c=.3",xlab="",ylab="",centers=TRUE)
#' set.seed(1)
#' plotCSarcs1D(Xp,Yp,t=1.5,c=.5,jit,main="t=1.5, c=.5",xlab="",ylab="",centers=TRUE)
#' set.seed(1)
#' plotCSarcs1D(Xp,Yp,t=2,c=.5,jit,main="t=2, c=.5",xlab="",ylab="",centers=TRUE)
#'
#' @export plotCSarcs1D
plotCSarcs1D <- function(Xp,Yp,t,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,centers=FALSE, ...)
{
  arcs.mid<-arcsCSmid.int(Xp,Yp,t,c)
  arcs.end<-arcsCSend.int(Xp,Yp,t)
  S<-c(arcs.mid$S, arcs.end$S)
  E<-c(arcs.mid$E, arcs.end$E)

  if (is.null(xlim))
  {xlim<-range(Xp,Yp)}

  ns<-length(S)
  jit<-(xlim[2]-xlim[1])*Jit
  ifelse(ns<=1,yjit<-0, yjit<-runif(ns,-jit,jit))

  if (is.null(ylim))
  {ylim<-2*c(-jit,jit)}

  if (is.null(main))
  {
    main.text=paste("Arcs of CS-PCD with t = ",t," and c = ",c,sep="")
    if (!centers){
      ifelse(ns<=1,main<-main.text,main<-c(main.text,"\n (arcs jittered along y-axis)"))
    } else
    {
      ifelse(ns<=1,main<-c(main.text,"\n (center added)"),main<-c(main.text,"\n (arcs jittered along y-axis & center added)"))
    }
  }

  nx<-length(Xp)
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
#NCS Functions that work for both middle and end intervals
#################################################################

#' @title The end points of the Central Similarity (CS) Proximity Region for a point - one interval case
#'
#' @description Returns the end points of the interval which constitutes the CS proximity region for a point in the
#' interval \code{int}\eqn{=(a,b)=}\code{(rv=1,rv=2)}.
#'
#' CS proximity region is constructed with respect to the interval \code{int}
#' with expansion parameter \eqn{t>0} and centrality parameter \eqn{c \in (0,1)}.
#'
#' Vertex regions are based on the (parameterized) center, \eqn{M_c},
#' which is \eqn{M_c=a+c(b-a)} for the interval, \code{int}\eqn{=(a,b)}.
#' The CS proximity region is constructed whether \code{x} is inside or outside the interval \code{int}.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param x A 1D point for which CS proximity region is constructed.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return The interval which constitutes the CS proximity region for the point \code{x}
#'
#' @seealso \code{\link{NPEint}} and \code{\link{NCStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' t<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' NCSint(7,int,t,c)
#' NCSint(17,int,t,c)
#' NCSint(1,int,t,c)
#' NCSint(-1,int,t,c)
#'
#' NCSint(3,int,t,c)
#' NCSint(4,int,t,c)
#' NCSint(a,int,t,c)
#'
#' @export NCSint
NCSint <- function(x,int,t,c=.5)
{
  if (!is.point(x,1) )
  {stop('x must be a scalar')}

  if (!is.point(t,1) || t <=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>y2)
  {stop('interval is degenerate or void, left end must be smaller than or equal to right end')}

  if (x<y1 || x>y2)
  {
    ifelse(x<y1,reg<-c(x-t*(y1-x),min(y1,x+t*(y1-x)) ),reg<-c( max(y2,x-t*(x-y2)),x+t*(x-y2)))
  } else
  {
    Mc<-y1+c*(y2-y1)
    if (x<=Mc)
    {
      reg <-c(max(y1,x-t*(x-y1)),min(y2,x+t*(1-c)*(x-y1)/c) )
    } else
    {
      reg<-c(max(y1,x-c*t*(y2-x)/(1-c)),min(y2,x+t*(y2-x)) )
    }
  }
  reg #proximity region interval
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for
#' Central Similarity  Proximity Catch Digraphs (CS-PCDs) - one interval case
#'
#' @description Returns \eqn{I(p_2} in \eqn{N_{CS}(p_1,t,c))} for points \eqn{p_1} and \eqn{p_2}, that is, returns 1 if \eqn{p_2} is in \eqn{N_{CS}(p_1,t,c)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t,c)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{t>0}
#' and centrality parameter \eqn{c \in (0,1)}.
#'
#' CS proximity region is constructed with respect to the
#' interval \eqn{(a,b)}. This function works whether \eqn{p_1} and \eqn{p_2} are inside or outside the interval \code{int}.
#'
#' Vertex regions for middle intervals are based on the center associated with the centrality parameter
#' \eqn{c \in (0,1)}. If \eqn{p_1} and \eqn{p_2} are identical, then it returns 1 regardless of their locations
#' (i.e., loops are allowed in the digraph).
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param p1 A 1D point for which the proximity region is constructed.
#' @param p2 A 1D point for which it is checked whether it resides in the proximity region
#' of \eqn{p_1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return \eqn{I(p_2} in \eqn{N_{CS}(p_1,t,c))} for p2, that is, returns 1 if \eqn{p_2} in \eqn{N_{CS}(p_1,t,c)}, returns 0 otherwise
#'
#' @seealso \code{\link{IarcCSmid.int}}, \code{\link{IarcCSend.int}} and \code{\link{IarcPEint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' t<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' IarcCSint(7,5,int,t,c)
#' IarcCSint(17,17,int,t,c)
#' IarcCSint(15,17,int,t,c)
#' IarcCSint(1,3,int,t,c)
#'
#' IarcCSint(-17,17,int,t,c)
#'
#' IarcCSint(3,5,int,t,c)
#' IarcCSint(3,3,int,t,c)
#' IarcCSint(4,5,int,t,c)
#' IarcCSint(a,5,int,t,c)
#'
#' c<-.4
#' r<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' IarcCSint(7,5,int,t,c)
#'
#' @export IarcCSint
IarcCSint <- function(p1,p2,int,t,c=.5)
{
  if (!is.point(p2,1) )
  {stop('p2 must be a scalar')}

  arc<-0
  pr<-NCSint(p1,int,t,c)  #proximity region as interval
  if (p2>=pr[1] && p2<=pr[2])
  {arc<-1}
  arc
} #end of the function
#'

#################################################################

#' @title The plot of the Central Similarity (CS) Proximity Regions for a general interval
#' (vertices jittered along \eqn{y}-coordinate) - one interval case
#'
#' @description Plots the points in and outside of the interval \code{int} and also the CS proximity regions (which are also intervals).
#' CS proximity regions are constructed with expansion parameter \eqn{t > 0} and centrality parameter \eqn{c \in (0,1)}.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default is \eqn{Jit=.1}) times range of proximity
#' regions and \code{Xp}) is added to the \eqn{y}-direction.
#' \code{center} is a logical argument, if \code{TRUE}, plot includes the
#' center of the interval as a vertical line in the plot, else center of the interval is not plotted.
#'
#' @param Xp A set of 1D points for which CS proximity regions are to be constructed.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis where \code{Jit} equals to the range of \code{Xp} and proximity region
#' intervals multiplied by \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges.
#' @param center A logical argument, if \code{TRUE}, plot includes the center of the interval
#' as a vertical line in the plot, else center of the interval is not plotted.
#' @param \dots	Additional \code{plot} parameters.
#'
#' @return Plot of the CS proximity regions for 1D points in or outside the interval \code{int}
#'
#' @seealso \code{\link{plotCSregs1D}}, \code{\link{plotCSregs}}, and \code{\link{plotPEregs.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' tau<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-10
#' xf<-(int[2]-int[1])*.1
#'
#' Xp<-runif(n,a-xf,b+xf)  #try also Xp<-runif(n,a-5,b+5)
#'
#' plotCSregs.int(7,int,tau,c,xlab="x",ylab="")
#'
#' plotCSregs.int(Xp,int,tau,c,xlab="x",ylab="")
#'
#' plotCSregs.int(17,int,tau,c,xlab="x",ylab="")
#' plotCSregs.int(1,int,tau,c,xlab="x",ylab="")
#' plotCSregs.int(4,int,tau,c,xlab="x",ylab="")
#'
#' plotCSregs.int(-7,int,tau,c,xlab="x",ylab="")
#'
#' @export plotCSregs.int
plotCSregs.int <- function(Xp,int,t,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,center=FALSE, ...)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}
  n<-length(Xp)
  pr<-c()
  for (i in 1:n)
  {x1<-Xp[i]
  pr<-rbind(pr,NCSint(x1,int,t,c))
  }

  if (is.null(xlim))
  {xlim<-range(Xp,int,pr)}

  xr<-xlim[2]-xlim[1]
  jit<-xr*Jit
  ifelse(n==1,yjit<-rep(0,n),yjit<-runif(n,-jit,jit))

  if (is.null(ylim))
  {ylim<-2*c(-jit,jit)}

  if (is.null(main))
  {
    main.text=paste("CS Proximity Regions with t = ",t," and c = ",c,sep="")
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

#' @title Number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' and quantities related to the interval - one interval case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' whose vertices are the
#' data points in \code{Xp} in the one middle interval case.
#' It also provides number of vertices
#' (i.e., number of data points inside the intervals)
#' and indices of the data points that reside in the intervals.
#'
#' The data points could be inside or outside the interval is \code{int}\eqn{=(a,b)}.
#'
#' CS proximity region is constructed
#' with an expansion parameter \eqn{t>0} and a centrality parameter \eqn{c \in (0,1)}.
#' CS proximity region is constructed for both points inside and outside the interval, hence
#' the arcs may exist for all points inside or outside the interval.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set of 1D points which constitute the vertices of CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and quantities related to the interval}
#' \item{num.arcs}{Total number of arcs in all intervals (including the end intervals),
#' i.e., the number of arcs for the entire CS-PCD}
#' \item{num.in.range}{Number of \code{Xp} points in the interval \code{int}}
#' \item{num.in.ints}{The vector of number of \code{Xp} points in the partition intervals (including the end intervals)}
#' \item{int.num.arcs}{The \code{vector} of the number of arcs of the component of the CS-PCD in the
#' partition intervals (including the end intervals)}
#' \item{data.int.ind}{A \code{vector} of indices of partition intervals in which data points reside.
#' Partition intervals are numbered from left to right with 1 being the left end interval.}
#' \item{ind.left.end, ind.mid, ind.right.end}{Indices of data points in the left end interval,
#' middle interval, and right end interval (respectively)}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support interval.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsCSmid.int}}, \code{\link{num.arcsCSend.int}},
#' and \code{\link{num.arcsPEint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' t<-2
#' a<-0; b<-10; int<-c(a,b)
#'
#' n<-10
#' set.seed(1)
#' Xp<-runif(n,a,b)
#' Narcs = num.arcsCSint(Xp,int,t,c)
#' Narcs
#' summary(Narcs)
#' plot(Narcs)
#'
#' @export num.arcsCSint
num.arcsCSint <- function(Xp,int,t,c=.5)
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
    ind.mid = which(Xp>=y1 & Xp <= y2)
    dat.mid<-Xp[ind.mid] # Xp points inside the int
    dat.left= Xp[Xp<y1]; dat.right= Xp[Xp>y2]
    ind.left.end = which(Xp<y1)
    ind.right.end = which(Xp>y2)
    int.ind[ind.left.end]=1
    int.ind[ind.mid]=2
    int.ind[ind.right.end]=3

    # number of arcs for the intervals
    narcs.left = num.arcsCSend.int(dat.left,int,t)
    narcs.right = num.arcsCSend.int(dat.right,int,t)
    narcs.mid = num.arcsCSmid.int(dat.mid,int,t,c)
    arcs = c(narcs.left,narcs.mid,narcs.right)

    ni.vec = c(length(dat.left),length(dat.mid),length(dat.right))

    narcs = sum(arcs)
  }

  NinInt = ni.vec[2]

  desc<-"Number of Arcs of the CS-PCD with vertices Xp and Quantities Related to the Support Interval"

  res<-list(desc=desc, #description of the output
            num.arcs=narcs, #number of arcs for the CS-PCD
            int.num.arcs=arcs, #vector of number of arcs for the partition intervals
            num.in.range=NinInt, #number of Xp points in the interval, int
            num.in.ints=ni.vec, #vector of numbers of Xp points in the partition intervals
            data.int.ind=int.ind, #indices of partition intervals in which data points reside, i.e., column number of part.int for each Xp point
            ind.mid =ind.mid, #indices of data points in the middle interval
            ind.left.end =ind.left.end, #indices of data points in the left end interval
            ind.right.end =ind.right.end, #indices of data points in the right end interval
            tess.points=int, #tessellation points
            vertices=Xp #vertices of the digraph
  )

  class(res)<-"NumArcs"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' and related quantities of the induced subdigraphs for points in the partition intervals -
#' multiple interval case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs and various other quantities related to the partition intervals
#' for Central Similarity Proximity Catch Digraph
#' (CS-PCD) whose vertices are the data points in \code{Xp}
#' in the multiple interval case.
#'
#' For this function, CS proximity regions are constructed data points inside or outside the intervals based
#' on \code{Yp} points with expansion parameter \eqn{t \ge 0} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points in the middle or end intervals.
#'
#' Range (or convex hull) of \code{Yp} (i.e., the interval \eqn{(\min(Yp),\max(Yp))}) is partitioned by the spacings based on
#' \code{Yp} points (i.e., multiple intervals are these partition intervals based on the order statistics of \code{Yp} points
#' whose union constitutes the range of \code{Yp} points). For the number of arcs, loops are not counted.
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the CS-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the partition intervals.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region;
#' must be \eqn{> 0}.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside the middle (partition) intervals
#' with the default \code{c=.5}.
#' For an interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and related quantities for the induced subdigraphs in the partition intervals}
#' \item{num.arcs}{Total number of arcs in all intervals (including the end intervals),
#' i.e., the number of arcs for the entire CS-PCD}
#' \item{num.in.range}{Number of \code{Xp} points in the range or convex hull of \code{Yp} points}
#' \item{num.in.ints}{The vector of number of \code{Xp} points in the partition intervals (including the end intervals)
#' based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the lengths of the middle partition intervals (i.e., end intervals excluded)
#' based on \code{Yp} points}
#' \item{int.num.arcs}{The \code{vector} of the number of arcs of the component of the CS-PCD in the
#' partition intervals (including the end intervals) based on \code{Yp} points}
#' \item{part.int}{A list of partition intervals based on \code{Yp} points}
#' \item{data.int.ind}{A \code{vector} of indices of partition intervals in which data points reside,
#' i.e., column number of \code{part.int} is provided for each \code{Xp} point. Partition intervals are numbered from left to right
#' with 1 being the left end interval.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the partition intervals based on \code{Yp} points.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsCSint}}, \code{\link{num.arcsCSmid.int}},
#' \code{\link{num.arcsCSend.int}}, and \code{\link{num.arcsPE1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' tau<-1.5
#' c<-.4
#' a<-0; b<-10; int<-c(a,b);
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
#' Narcs = num.arcsCS1D(Xp,Yp,tau,c)
#' Narcs
#' summary(Narcs)
#' plot(Narcs)
#'
#' @export num.arcsCS1D
num.arcsCS1D <- function(Xp,Yp,t,c=.5)
{
  if (!is.point(Yp,length(Yp)) )
  {stop('Yp must be 1D vector of numerical entries')}

  nx<-length(Xp); ny<-length(Yp)

  if (nx==0 || ny==0)
  {stop('Not enough points to construct CS-PCD')}

  Ys<-sort(Yp)  #sorted data points from classes X and Y
  ymin<-Ys[1]; ymax<-Ys[ny];
  Yrange=c(ymin, ymax)

  int.ind = rep(0,nx)
  dat.mid<-Xp[Xp>ymin & Xp < ymax] # Xp points inside  min(Yp) and max (Yp)
  dat.left= Xp[Xp<ymin]; dat.right= Xp[Xp>ymax]
  int.ind[which(Xp<ymin)]=1
  int.ind[which(Xp>ymax)]=ny+1

  #for end intervals
  narcs.left = num.arcsCSend.int(dat.left,Yrange,t)
  narcs.right = num.arcsCSend.int(dat.right,Yrange,t)
  arcs=narcs.left

  #for middle intervals
  n.int<-ny-1 #number of Yp middle intervals
  nx2<-length(dat.mid)  #number of Xp points inside the middle intervals

  Wvec<-Yspacings<-vector()
  for (j in 1:n.int)
  {
    Yspacings = rbind(Yspacings,c(Ys[j],Ys[j+1]))
    Wvec<-c(Wvec,Ys[j+1]-Ys[j])
  }

  part.ints =rbind(c(-Inf,ymin),Yspacings,c(ymax,Inf))

  ni.vec = vector()
  for (i in 1:n.int)
  {
    ind = which(Xp>=Ys[i] & Xp < Ys[i+1])
    dat.int<-Xp[ind] #X points in the ith Yp mid interval
    int.ind[ind] = i+1
    ni.vec = c(ni.vec,length(dat.int))
    narcs.mid = num.arcsCSmid.int(dat.int,Yspacings[i,],t,c)
    arcs = c(arcs,narcs.mid)
  }
  ni.vec = c(length(dat.left),ni.vec,length(dat.right))
  arcs = c(arcs,narcs.right) #reordering the number of arcs vector according to the order of the intervals, from left to right
  narcs = sum(arcs)

  desc<-"Number of Arcs of the CS-PCD with vertices Xp and Related Quantities for the Induced Subdigraphs for the Points in the Partition Intervals"

  res<-list(desc=desc, #description of the output
            num.arcs=narcs, #number of arcs for the entire PCD
            int.num.arcs=arcs, #vector of number of arcs for the partition intervals
            num.in.range=nx2, #number of Xp points in the range of Yp points
            num.in.ints=ni.vec, #number of Xp points in the partition intervals
            weight.vec=Wvec, #lengths of the middle partition intervals
            partition.intervals=t(part.ints), #matrix of the partition intervals, each column is one interval
            data.int.ind=int.ind, #indices of partition intervals in which data points reside, i.e., column number of part.int for each Xp point
            tess.points=Yp, #tessellation points
            vertices=Xp #vertices of the digraph
  )

  class(res)<-"NumArcs"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title The arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for 1D data - one interval case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of CS-PCD.  \code{int} determines the end points of the interval.
#'
#' For this function, CS proximity regions are constructed data points inside or outside the interval based
#' on \code{int} points with expansion parameter \eqn{t > 0} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points in the middle or end intervals.
#' It also provides various descriptions and quantities about the arcs of the CS-PCD
#' such as number of arcs, arc density, etc.
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the CS-PCD.
#' @param int A \code{vector} of two 1D points which constitutes the end points of the interval.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
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
#' \item{S}{Tails (or sources) of the arcs of CS-PCD for 1D data}
#' \item{E}{Heads (or arrow ends) of the arcs of CS-PCD for 1D data}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsCS1D}}, \code{\link{arcsCSmid.int}}, \code{\link{arcsCSend.int}}, and \code{\link{arcsPE1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' tau<-2
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
#' Arcs<-arcsCSint(Xp,int,tau,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' Xp<-runif(n,a+10,b+10)
#' Arcs=arcsCSint(Xp,int,tau,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' @export arcsCSint
arcsCSint <- function(Xp,int,t,c=.5)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(int))

  if (!is.point(Xp,length(Xp)) || !is.point(int,length(int)) )
  {stop('Xp and int must be 1D vectors of numerical entries')}

  if (!is.point(t,1) || t <=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp)
  S<-E<-vector()  #S is for source and E is for end points for the arcs
  if (nx==0)
  {stop('Not enough points to construct CS-PCD')}

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
      xR<-x1+t*(y1-x1); xL<-x1-t*(y1-x1)
      ind.tails<-((xLe < min(xR,y1)) & (xLe > xL))
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
      v<-rel.vert.mid.int(x1,int,c)$rv
      if (v==1)
      {
        xR<-x1+t*(1-c)*(x1-y1)/c
        xL<-x1-t*(x1-y1)
        ind.tails<-((Xinl < min(xR,y2)) & (Xinl > max(xL,y1)))
        st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
        S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
      } else {
        xR <-x1+t*(y2-x1)
        xL <-x1-c*t*(y2-x1)/(1-c)
        ind.tails<-((Xinl < min(xR,y2)) & (Xinl > max(xL,y1)))
        st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
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
      xR<-x1+t*(x1-y2); xL<-x1-t*(x1-y2)
      ind.tails<-((xRe < xR) & xRe > max(y2,xL))
      st<-sum(ind.tails)  #sum of tails of the arcs with head XRe[j]
      S<-c(S,rep(x1,st)); E<-c(E,xRe[ind.tails])
      }
    }
  }

  if (length(S)==0)
  {S<-E<-NA}

  param<-list(c,t)
  names(param)<-c("centrality parameter","expansion parameter")

  typ<-paste("Central Similarity Proximity Catch Digraph (CS-PCD) for 1D Points with Expansion Parameter t = ",t, " and Centrality Parameter c = ",c,sep="")

  main.txt<-paste("Arcs of CS-PCD with t = ",round(t,2)," and c = ",round(c,2),"\n (arcs jittered along y-axis)",sep="")

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

#' @title The plot of the arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs) for 1D data
#' (vertices jittered along \eqn{y}-coordinate) - one interval case
#'
#' @description Plots the arcs of CS-PCD whose vertices are the 1D points, \code{Xp}. CS proximity regions are constructed with
#' expansion parameter \eqn{t > 0} and centrality parameter \eqn{c \in (0,1)} and the intervals are based on the
#' interval \code{int}\eqn{=(a,b)}
#' That is, data set \code{Xp}
#' constitutes the vertices of the digraph and \code{int} determines the end points of the interval.
#'
#' For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)} (default for \eqn{Jit=.1}) is added to
#' the \eqn{y}-direction where \code{Jit} equals to the range of \eqn{\{}\code{Xp}, \code{int}\eqn{\}}
#' multiplied by \code{Jit} with default for \eqn{Jit=.1}).
#' \code{center} is a logical argument, if \code{TRUE}, plot includes the center of the interval \code{int}
#' as a vertical line in the plot, else center of the interval is not plotted.
#'
#' @param Xp A \code{vector} of 1D points constituting the vertices of the CS-PCD.
#' @param int A \code{vector} of two 1D points constituting the end points of the interval.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
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
#' @return A plot of the arcs of CS-PCD whose vertices are the 1D data set \code{Xp} in which vertices are jittered
#' along \eqn{y}-axis for better visualization.
#'
#' @seealso \code{\link{plotCSarcs1D}} and \code{\link{plotPEarcs.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' tau<-2
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
#' Ylim=3*c(-1,1)
#'
#' jit<-.1
#' plotCSarcs.int(Xp,int,t=tau,c,jit,xlab="",ylab="",xlim=Xlim,ylim=Ylim)
#'
#' set.seed(1)
#' plotCSarcs.int(Xp,int,t=1.5,c=.3,jit,xlab="",ylab="",center=TRUE)
#' set.seed(1)
#' plotCSarcs.int(Xp,int,t=2,c=.4,jit,xlab="",ylab="",center=TRUE)
#'
#' @export plotCSarcs.int
plotCSarcs.int <- function(Xp,int,t,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,center=FALSE, ...)
{
  arcs<-arcsCS1D(Xp,int,t,c)
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
    main.text=paste("Arcs of CS-PCD with t = ",t," and c = ",c,sep="")
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

#' @title The arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for 1D data - multiple interval case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for 1D data set \code{Xp} as the vertices
#' of CS-PCD and related parameters and the quantities of the digraph.
#' \code{Yp} determines the end points of the intervals.
#'
#' For this function, CS proximity regions are constructed data points inside or outside the intervals based
#' on \code{Yp} points with expansion parameter \eqn{t>0} and centrality parameter \eqn{c \in (0,1)}. That is, for this function,
#' arcs may exist for points in the middle or end intervals.
#' It also provides various descriptions and quantities about the arcs of the CS-PCD
#' such as number of arcs, arc density, etc.
#'
#' Equivalent to function \code{\link{arcsCS1D}}.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points which constitute the vertices of the CS-PCD.
#' @param Yp A set or \code{vector} of 1D points which constitute the end points of the intervals.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, they are expansion and centrality parameters.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the intervalization of the real line based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of CS-PCD for 1D data}
#' \item{E}{Heads (or arrow ends) of the arcs of CS-PCD for 1D data}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsCSend.int}}, \code{\link{arcsCSmid.int}}, \code{\link{arcsCS1D}}, and \code{\link{arcsPE1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-2
#' c<-.4
#' a<-0; b<-10;
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xr<-range(a,b)
#' xf<-(xr[2]-xr[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' Arcs<-arcsCS1D(Xp,Yp,t,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' S<-Arcs$S
#' E<-Arcs$E
#'
#' arcsCS1D(Xp,Yp,t,c)
#'
#' arcsCS1D(Xp,Yp+10,t,c)
#'
#' jit<-.1
#' yjit<-runif(nx,-jit,jit)
#'
#' Xlim<-range(a,b,Xp,Yp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),
#' main="arcs of CS-PCD for points (jittered along y-axis)\n in middle intervals ",
#' xlab=" ", ylab=" ", xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit),pch=".")
#' abline(h=0,lty=1)
#' points(Xp, yjit,pch=".",cex=3)
#' abline(v=Yp,lty=2)
#' arrows(S, yjit, E, yjit, length = .05, col= 4)
#'
#' t<-2
#' c<-.4
#' a<-0; b<-10;
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' arcsCS1D(Xp,Yp,t,c)
#'
#' @export arcsCS1D
arcsCS1D <- function(Xp,Yp,t,c=.5)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)
  S<-E<-vector()  #S is for source and E is for end points for the arcs
  if (nx==0 || ny==0)
  {stop('Not enough points to construct CS-PCD')}

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
      xR<-x1+t*(ymin-x1); xL<-x1-t*(ymin-x1)
      ind.tails<-((xLe < min(xR,ymin)) & (xLe > xL))
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
          v<-rel.vert.mid.int(x1,int,c)$rv
          if (v==1)
          {
            xR<-x1+t*(1-c)*(x1-y1)/c
            xL<-x1-t*(x1-y1)
            ind.tails<-((Xinl < min(xR,y2)) & (Xinl > max(xL,y1)))
            st<-sum(ind.tails)  #sum of tails of the arcs with head Xi[j]
            S<-c(S,rep(x1,st)); E<-c(E,Xinl[ind.tails])
          } else {
            xR <-x1+t*(y2-x1)
            xL <-x1-c*t*(y2-x1)/(1-c)
            ind.tails<-((Xinl < min(xR,y2)) & (Xinl > max(xL,y1)))
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
      xR<-x1+t*(x1-ymax); xL<-x1-t*(x1-ymax)
      ind.tails<-((xRe < xR) & xRe > max(ymax,xL))
      st<-sum(ind.tails)  #sum of tails of the arcs with head XRe[j]
      S<-c(S,rep(x1,st)); E<-c(E,xRe[ind.tails])
      }
    }
  }
  if (length(S)==0)
  {S<-E<-NA}

  param<-c(c,t)
  names(param)<-c("centrality parameter","expansion parameter")

  typ<-paste("Central Similarity Proximity Catch Digraph (CS-PCD) for 1D Points with Expansion Parameter t = ",t, " and Centrality Parameter c = ",c,sep="")

  main.txt<-paste("Arcs of CS-PCD with t = ",round(t,2)," and c = ",round(c,2),"\n (arcs jittered along y-axis)",sep="")


  nvert<-nx; nint<-ny+1; narcs=ifelse(!is.null(S),length(S),0);
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

#' @title Incidence matrix for Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' for 1D data - multiple interval case
#'
#' @description Returns the incidence matrix for the CS-PCD for a given 1D numerical data set, \code{Xp},
#' as the vertices of the digraph and \code{Yp} determines the end points of the intervals (in the multi-interval case).
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' CS proximity region is constructed
#' with an expansion parameter \eqn{t>0} and a centrality parameter \eqn{c \in (0,1)}.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp a set of 1D points which constitutes the vertices of the digraph.
#' @param Yp a set of 1D points which constitutes the end points of the intervals
#' that partition the real line.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return Incidence matrix for the CS-PCD with vertices being 1D data set, \code{Xp},
#' and \code{Yp} determines the end points of the intervals (the multi-interval case)
#'
#' @seealso \code{\link{inci.matCS1D}}, \code{\link{inci.matPEtri}}, and \code{\link{inci.matPE}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-2
#' c<-.4
#' a<-0; b<-10;
#' nx<-10; ny<-4
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' IM<-inci.matCS1D(Xp,Yp,t,c)
#' IM
#' dom.num.greedy(IM)
#' \donttest{
#' dom.num.exact(IM)  #might take a long time depending on nx
#' }
#' Idom.num.up.bnd(IM,5)
#'
#' Arcs<-arcsCS1D(Xp,Yp,t,c)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' inci.matCS1D(Xp,Yp+10,t,c)
#'
#' t<-2
#' c<-.4
#' a<-0; b<-10;
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' inci.matCS1D(Xp,Yp,t,c)
#'
#' @export inci.matCS1D
inci.matCS1D <- function(Xp,Yp,t,c=.5)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)
  nt<-ny-1 #number of Yp middle intervals
  if (nx==0 || ny==0)
  {stop('Not enough points to construct CS-PCD')}

  if (nx>=1)
  {
    Ys<-sort(Yp)  #sorted data points from classes X and Y
    ymin<-Ys[1]; ymax<-Ys[ny];

    pr<-c()
    for (i in 1:nx)
    { x1<-Xp[i]
    if (x1<ymin || x1>=ymax)
    {int<-c(ymin,ymax)
    pr<-rbind(pr,NCSint(x1,int,t,c))
    }
    if (nt>=1)
    {
      for (j in 1:nt)
      {
        if (x1>=Ys[j] & x1 < Ys[j+1] )
        {  y1<-Ys[j]; y2<-Ys[j+1]; int<-c(y1,y2)
        pr<-rbind(pr,NCSint(x1,int,t,c))
        }
      }
    }
    }

    inc.mat<-matrix(0, nrow=nx, ncol=nx)
    for (i in 1:nx)
    { reg<-pr[i,]
    for (j in 1:nx)
    {
      inc.mat[i,j]<-sum(Xp[j]>=reg[1] & Xp[j]<reg[2])
    }
    }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title Incidence matrix for Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' for 1D data - one interval case
#'
#' @description Returns the incidence matrix for the CS-PCD for a given 1D numerical data set, \code{Xp},
#' as the vertices of the digraph and \code{int} determines the end points of the interval (in the one interval case).
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' CS proximity region is constructed
#' with an expansion parameter \eqn{t > 0} and a centrality parameter \eqn{c \in (0,1)}.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp a set of 1D points which constitutes the vertices of the digraph.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return Incidence matrix for the CS-PCD with vertices being 1D data set, \code{Xp},
#' and \code{int} determines the end points of the intervals (in the one interval case)
#'
#' @seealso \code{\link{inci.matCS1D}}, \code{\link{inci.matPE1D}},
#' \code{\link{inci.matPEtri}}, and \code{\link{inci.matPE}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' t<-1
#' a<-0; b<-10; int<-c(a,b)
#'
#' xf<-(int[2]-int[1])*.1
#'
#' set.seed(123)
#'
#' n<-10
#' Xp<-runif(n,a-xf,b+xf)
#'
#' IM<-inci.matCSint(Xp,int,t,c)
#' IM
#'
#' dom.num.greedy(IM)
#' Idom.num.up.bnd(IM,3)
#' dom.num.exact(IM)
#'
#' inci.matCSint(Xp,int+10,t,c)
#' }
#'
#' @export inci.matCSint
inci.matCSint <- function(Xp,int,t,c=.5)
{
  if (!is.point(Xp,length(Xp)) )
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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
    pr<-rbind(pr,NCSint(x1,int,t,c))
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

#' @title The plot of the Central Similarity (CS) Proximity Regions (vertices jittered along \eqn{y}-coordinate)
#' - multiple interval case
#'
#' @description Plots the points in and outside of the intervals based on \code{Yp} points and also the CS proximity regions
#' (which are also intervals).
#'
#' CS proximity region is constructed with expansion parameter \eqn{t>0} and
#' centrality parameter \eqn{c \in (0,1)}. For better visualization, a uniform jitter from \eqn{U(-Jit,Jit)}
#' (default is \eqn{Jit=.1}) times range of \code{Xp} and \code{Yp} and the proximity regions (intervals)) is added to the
#' \eqn{y}-direction.
#'
#' \code{centers} is a logical argument, if \code{TRUE},
#' plot includes the centers of the intervals as vertical lines in the plot,
#' else centers of the intervals are not plotted.
#'
#' See also (\insertCite{ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Xp A set of 1D points for which CS proximity regions are plotted.
#' @param Yp A set of 1D points which constitute the end points of the intervals which
#' partition the real line.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside middle intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param Jit A positive real number that determines the amount of jitter along the \eqn{y}-axis, default=\code{0.1} and
#' \code{Xp} points are jittered according to \eqn{U(-Jit,Jit)} distribution along the \eqn{y}-axis where \code{Jit} equals to the range of \code{Xp} and \code{Yp} and the
#' proximity regions (intervals) multiplied by \code{Jit}).
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles of the \eqn{x} and \eqn{y} axes in the plot (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param centers A logical argument, if \code{TRUE}, plot includes the centers of the intervals
#' as vertical lines in the plot, else centers of the intervals are not plotted.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the CS proximity regions for 1D points located in the middle or end intervals
#' based on \code{Yp} points
#'
#' @seealso \code{\link{plotCSregs.int}} and \code{\link{plotPEregs1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-2
#' c<-.4
#' a<-0; b<-10;
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' xr<-range(a,b)
#' xf<-(xr[2]-xr[1])*.1
#'
#' Xp<-runif(nx,a-xf,b+xf)
#' Yp<-runif(ny,a,b)
#'
#' plotCSregs1D(Xp,Yp,t,c,xlab="",ylab="")
#'
#' plotCSregs1D(Xp,Yp+10,t,c,xlab="",ylab="")
#'
#' @export plotCSregs1D
plotCSregs1D <- function(Xp,Yp,t,c=.5,Jit=.1,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,centers=FALSE, ...)
{
  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
  {stop('c must be a scalar in (0,1)')}

  nx<-length(Xp); ny<-length(Yp)

  if (ny < 1 || nx < 1)
  {stop('Both Xp and Yp points must be nonempty to construct PE-PCD')}

  LE<-RE<-vector()
  if (nx>=1)
  { Xp<-sort(Xp)
  Ys<-sort(Yp)  #sorted data points from classes X and Y
  ymin<-Ys[1]; ymax<-Ys[ny];

  ifelse(nx==1,yjit<-0,yjit<-runif(nx,-Jit,Jit))

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
    pr<-NCSint(x1,int,t,c)
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
        pr<-NCSint(x1,int,t,c)
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
    pr<-NCSint(x1,int,t,c)
    LE<-c(LE,pr[1]); RE<-c(RE,pr[2])
    }
  }
  }

  if (is.null(xlim))
  {xlim<-range(Xp,Yp,LE,RE)}
  if (is.null(ylim))
  {ylim<-2*c(-Jit,Jit)}

  xr<-xlim[2]-xlim[1]

  if (is.null(main))
  {
    main.text=paste("CS Proximity Regions with t = ",t," and c = ",c,sep="")
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

#' @title The indicator for a point being a dominating point for Central Similarity
#' Proximity Catch Digraphs (CS-PCDs) for an interval
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of CS-PCD\eqn{)} where the vertices of the CS-PCD are the 1D data set \code{Xp}).
#'
#' CS proximity region is defined with respect to the interval \code{int} with an expansion parameter, \eqn{t>0},
#' and a centrality parameter, \eqn{c \in (0,1)}, so arcs may exist for \code{Xp} points inside the interval \code{int}\eqn{=(a,b)}.
#'
#' Vertex regions are based on the center associated with the centrality parameter \eqn{c \in (0,1)}.
#' \code{rv} is the index of the vertex region \code{p} resides, with default=\code{NULL}.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' @param p A 1D point that is to be tested for being a dominating point or not of the CS-PCD.
#' @param Xp A set of 1D points which constitutes the vertices of the CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param rv Index of the vertex region in which the point resides, either \code{1,2} or \code{NULL}
#' (default is \code{NULL}).
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of CS-PCD\eqn{)} where the vertices of the CS-PCD are the 1D data set \code{Xp}),
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num1PEint}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' t<-2
#' c<-.4
#' a<-0; b<-10; int<-c(a,b)
#'
#' Mc<-centerMc(int,c)
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif(n,a,b)
#'
#' Idom.num1CSint(Xp[5],Xp,int,t,c)
#'
#' Idom.num1CSint(2,Xp,int,t,c,ch.data.pnt = FALSE)
#' #gives an error if ch.data.pnt = TRUE since p is not a data point in Xp
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1CSint(Xp[i],Xp,int,t,c))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' domset<-Xp[ind.gam1]
#' if (length(ind.gam1)==0)
#' {domset<-NA}
#'
#' #or try
#' Rv<-rel.vert.mid.int(Xp[5],int,c)$rv
#' Idom.num1CSint(Xp[5],Xp,int,t,c,Rv)
#'
#' Xlim<-range(a,b,Xp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),xlab="",pch=".",xlim=Xlim+xd*c(-.05,.05))
#' abline(h=0)
#' abline(v=c(a,b,Mc),col=c(1,1,2),lty=2)
#' points(cbind(Xp,0))
#' points(cbind(domset,0),pch=4,col=2)
#' text(cbind(c(a,b,Mc),-0.1),c("a","b","Mc"))
#'
#' Idom.num1CSint(Xp[5],Xp,int,t,c)
#'
#' n<-10
#' Xp2<-runif(n,a+b,b+10)
#' Idom.num1CSint(5,Xp2,int,t,c)
#'
#' @export Idom.num1CSint
Idom.num1CSint <- function(p,Xp,int,t,c=.5,rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p,1) )
  {stop('p must be a scalar')}

  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(c,1) || c<=0 || c>=1)
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
  {rv<-rel.vert.mid.int(p,int,c)$rv #determines the vertex region for 1D point p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  Xp<-Xp[(Xp>=y1 & Xp<=y2)]
  n<-length(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (IarcCSint(p,Xp[i],int,t,c)==0)
  {dom<-0;}
    i<-i+1;
  }
  dom
} #end of the function
#'
