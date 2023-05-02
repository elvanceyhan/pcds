#AuxGeometry.R;
#Contains the ancillary functions for geometric constructs such as equation of lines for two points
#################################################################

#' @title Check the argument is a point of a given dimension
#'
#' @description Returns \code{TRUE} if the argument \code{p} is a \code{numeric} point of dimension \code{dim}
#' (default is \code{dim=2}); otherwise returns \code{FALSE}.
#'
#' @param p A \code{vector} to be checked to see it is a point of dimension \code{dim} or not.
#' @param dim A positive integer representing the dimension of the argument \code{p}.
#'
#' @return \code{TRUE} if \code{p} is a \code{vector} of dimension \code{dim}.
#'
#' @seealso \code{\link{dimension}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(-1.22,-2.33); B<-c(2.55,3.75,4)
#' is.point(A)
#' is.point(A,1)
#'
#' is.point(B)
#' is.point(B,3)
#' }
#'
#' @export
is.point <- function(p,dim=2)
{
  res<-is.numeric(p)==TRUE && is.vector(p)==TRUE && !is.list(p) && length(p)==dim
  res
}

#################################################################

#' @title The dimension of a \code{vector} or matrix or a data frame
#'
#' @description
#' Returns the dimension (i.e., number of columns) of \code{x}, which is a matrix or a \code{vector} or a data
#' frame. This is different than the \code{dim} function in the base distribution of \code{R}, in the sense that,
#' \code{dimension} gives only the number of columns of the argument \code{x}, while \code{dim} gives the number of rows and
#' columns of \code{x}.
#' \code{dimension} also works for a scalar or a vector, while \code{dim} yields \code{NULL} for such arguments.
#'
#' @param x A \code{vector} or a matrix or a data frame whose dimension is to be determined.
#'
#' @return Dimension (i.e., number of columns) of \code{x}
#'
#' @seealso \code{\link{is.point}} and \code{\link{dim}} from the base distribution of \code{R}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' dimension(3)
#' dim(3)
#'
#' A<-c(1,2)
#' dimension(A)
#' dim(A)
#'
#' B<-c(2,3)
#' dimension(rbind(A,B,A))
#' dimension(cbind(A,B,A))
#'
#' M<-matrix(runif(20),ncol=5)
#' dimension(M)
#' dim(M)
#'
#' dimension(c("a","b"))
#' }
#'
#' @export dimension
dimension <- function(x)
{
  if ( (!is.vector(x) && !is.matrix(x) && !is.data.frame(x) )  || (is.list(x) && !is.data.frame(x)) )
  {stop('x must be a vector or a matrix or a data frame')}

  if (is.vector(x))
  {dimen<-length(x)
  } else
  {x<-as.matrix(x)
  dimen<-ncol(x)
  }
  dimen
}


#################################################################

#' @title The slope of a line
#'
#' @description Returns the slope of the line joining two distinct 2D points \code{a} and \code{b}.
#'
#' @param a,b 2D points that determine the straight line (i.e., through which the straight line passes).
#'
#' @return Slope of the line joining 2D points \code{a} and \code{b}
#'
#' @seealso \code{\link{Line}}, \code{\link{paraline}}, and \code{\link{perpline}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' A<-c(-1.22,-2.33); B<-c(2.55,3.75)
#' slope(A,B)
#'
#' slope(c(1,2),c(2,3))
#'
#' @export slope
slope <- function(a,b)
{
  if (!is.point(a) || !is.point(b))
  {stop('both arguments must be numeric points (vectors) of dimension 2')}

  if (all(a==b))
  {stop('both arguments are same; slope is not well defined')}

  (b[2]-a[2])/(b[1]-a[1])
}

#################################################################

#' @title The line joining two distinct 2D points \code{a} and \code{b}
#'
#' @description
#' An object of class \code{"Lines"}.
#' Returns the \code{equation, slope, intercept}, and \eqn{y}-coordinates of the line crossing two
#' distinct 2D points \code{a} and \code{b} with \eqn{x}-coordinates provided in \code{vector} \code{x}.
#'
#' This function is different from the \code{\link[stats]{line}} function in the standard \code{stats} package
#' in \code{R} in the sense that \code{Line(a,b,x)} fits the line passing
#' through points \code{a} and \code{b} and returns various quantities (see below) for this line and \code{x} is
#' the \eqn{x}-coordinates of the points we want to find on the \code{Line(a,b,x)} while \code{line(a,b)} fits the line robustly
#' whose \eqn{x}-coordinates are in \code{a} and \eqn{y}-coordinates are in \code{b}.
#'
#' \code{Line(a,b,x)} and \code{line(x,Line(A,B,x)$y)}
#' would yield the same straight line (i.e., the line with the same coefficients.)
#'
#' @param a,b 2D points that determine the straight line (i.e., through which the straight line passes).
#' @param x A scalar or a \code{vector} of scalars representing the \eqn{x}-coordinates of the line.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A description of the line}
#' \item{mtitle}{The \code{"main"} title for the plot of the line}
#' \item{points}{The input points \code{a} and \code{b} through which the straight line passes
#' (stacked row-wise, i.e., row 1 is point \code{a} and row 2 is point \code{b}).}
#' \item{x}{The input scalar or \code{vector} which constitutes the \eqn{x}-coordinates of the point(s) of interest
#' on the line.}
#' \item{y}{The output scalar or \code{vector} which constitutes the \eqn{y}-coordinates of the point(s) of interest
#' on the line. If \code{x} is a scalar, then \code{y} will be a scalar and if \code{x} is a \code{vector} of scalars,
#' then \code{y} will be a \code{vector} of scalars.}
#' \item{slope}{Slope of the line, \code{Inf} is allowed, passing through points \code{a} and \code{b} }
#' \item{intercept}{Intercept of the line passing through points \code{a} and \code{b}}
#' \item{equation}{Equation of the line passing through points \code{a} and \code{b}}
#'
#' @seealso \code{\link{slope}}, \code{\link{paraline}}, \code{\link{perpline}}, \code{\link[stats]{line}}
#' in the generic \code{stats} package and and \code{\link{Line3D}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(-1.22,-2.33); B<-c(2.55,3.75)
#'
#' xr<-range(A,B);
#' xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#'
#' lnAB<-Line(A,B,x)
#' lnAB
#' summary(lnAB)
#' plot(lnAB)
#'
#' line(A,B)  #this takes vector A as the x points and vector B as the y points and fits the line
#' #for example, try
#' x=runif(100); y=x+(runif(100,-.05,.05))
#' plot(x,y)
#' line(x,y)
#'
#' x<-lnAB$x
#' y<-lnAB$y
#' Xlim<-range(x,A,B)
#' if (!is.na(y[1])) {Ylim<-range(y,A,B)} else {Ylim<-range(A,B)}
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' pf<-c(xd,-yd)*.025
#'
#' #plot of the line joining A and B
#' plot(rbind(A,B),pch=1,xlab="x",ylab="y",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' if (!is.na(y[1])) {lines(x,y,lty=1)} else {abline(v=A[1])}
#' text(rbind(A+pf,B+pf),c("A","B"))
#' int<-my_round(lnAB$intercep,2)  #intercept
#' sl<-my_round(lnAB$slope,2)  #slope
#' text(rbind((A+B)/2+pf*3),ifelse(is.na(int),paste("x=",A[1]),
#' ifelse(sl==0,paste("y=",int),
#' ifelse(sl==1,ifelse(sign(int)<0,paste("y=x",int),paste("y=x+",int)),
#' ifelse(sign(int)<0,paste("y=",sl,"x",int),paste("y=",sl,"x+",int))))))
#' }
#'
#' @export Line
Line <- function(a,b,x)
{
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))

  if (!is.point(a) || !is.point(b))
  {stop('a and b must be numeric points of dimension 2')}

  if (!is.point(x,length(x)))
  {stop('x must be a numeric vector')}

  if (all(a==b))
  {stop('a and b are same; line is not well defined')}

  sl<-slope(a,b)
  if (abs(sl)==Inf)
  {
    int<-ln<-NA
  } else
  {
    int<-(b[2]+(0-b[1])*sl);
    ln<-(b[2]+(x-b[1])*sl);
  }

  names(sl)<-"slope"
  names(int)<-"intercept"

  line.desc<-paste("Line Passing through two Distinct Points", aname, "and", bname)
  main.txt<-paste("Line Passing through Points", aname, "and", bname)

  pts<-rbind(a,b)
  row.names(pts)<-c(aname,bname)

  if (abs(sl)==Inf)
  { vert<-a[1]
  eqn<-reqn<-paste("x=",vert)
  } else
  {
    vert<-NA
    eqn<-ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),
                                            ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,
                                                                                                ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
    rsl<-my_round(sl,2)
    rint<-my_round(int,2)
    reqn<-ifelse(rsl==0 & rint==0,"y=0",ifelse(rsl!=0 & rint==0,paste("y=",rsl,"x",sep=""),
                                               ifelse(rsl==0 & rint!=0,paste("y=",rint,sep=""),ifelse(rsl==1,
                                                                                                      ifelse(sign(rint)<0,paste("y=x",rint,sep=""),paste("y=x+",rint,sep="")),
                                                                                                      ifelse(sign(rint)<0,paste("y=",rsl,"x",rint,sep=""),paste("y=",rsl,"x+",rint,sep="")))) ))
  }

  res<-list(
    desc=line.desc,
    mtitle=main.txt,
    points=pts,
    vert=vert,
    x=x,
    y=ln,
    slope=sl,
    intercept=int,
    equation=eqn,
    eqnlabel=reqn
  )
  class(res)<-"Lines"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The point of intersection of two lines defined by two pairs of points
#'
#' @description Returns the intersection of two lines, first line passing through points \code{p1} and \code{q1}
#' and second line passing through points \code{p2} and \code{q2}. The points are chosen so that lines are well
#' defined.
#'
#' @param p1,q1 2D points that determine the first straight line (i.e., through which the first straight line
#' passes).
#' @param p2,q2 2D points that determine the second straight line (i.e., through which the second straight line
#' passes).
#'
#' @return The coordinates of the point of intersection of the two lines, first passing through points
#' \code{p1} and \code{q1} and second passing through points \code{p2} and \code{q2}.
#'
#' @seealso \code{\link{int.circ.line}} and \code{\link{dist.point2line}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(-1.22,-2.33); B<-c(2.55,3.75); C<-c(0,6); D<-c(3,-2)
#'
#' ip<-int.2lines(A,B,C,D)
#' ip
#' pts<-rbind(A,B,C,D,ip)
#' xr<-range(pts[,1])
#' xf<-abs(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#' lnAB<-Line(A,B,x)
#' lnCD<-Line(C,D,x)
#'
#' y1<-lnAB$y
#' y2<-lnCD$y
#' Xlim<-range(x,pts)
#' Ylim<-range(y1,y2,pts)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' pf<-c(xd,-yd)*.025
#'
#' #plot of the line joining A and B
#' plot(rbind(A,B,C,D),pch=1,xlab="x",ylab="y",main="Point of Intersection of Two Lines",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' lines(x,y1,lty=1,col=1)
#' lines(x,y2,lty=1,col=2)
#' text(rbind(A+pf,B+pf),c("A","B"))
#' text(rbind(C+pf,D+pf),c("C","D"))
#' text(rbind(ip+pf),c("intersection\n point"))
#' }
#'
#' @export int.2lines
int.2lines <- function(p1,q1,p2,q2)
{
  if (!is.point(p1) || !is.point(p2) || !is.point(q1) || !is.point(q2))
  {stop('All arguments must be numeric points of dimension 2')}

  if (all(p1==q1) || all(p2==q2))
  {stop('At least one of the lines is not well defined')}

  sl1<-slope(p1,q1); sl2<-slope(p2,q2);
  if (sl1==sl2)
  {stop('The lines are parallel, hence do not intersect at a point')}

  p11<-p1[1]; p12<-p1[2]; p21<-p2[1]; p22<-p2[2];
  q11<-q1[1]; q12<-q1[2]; q21<-q2[1]; q22<-q2[2];
  if (p11==q11 && p12!=q12 && sl2!=Inf)
  {int.pnt<-c(p11,(p11*p22-p11*q22+p21*q22-p22*q21)/(-q21+p21))
  return(int.pnt);stop}

  if (p21==q21 && p22!=q22 && sl1!=Inf)
  {int.pnt<-c(p21,(p11*q12+p12*p21-p12*q11-p21*q12)/(-q11+p11))
  return(int.pnt);stop}

  int.pnt<-c((p11*p21*q12-p11*p21*q22+p11*p22*q21-p11*q12*q21-p12*p21*q11+p12*q11*q21+p21*q11*q22-p22*q11*q21)/(p11*p22-p11*q22-p12*p21+p12*q21+p21*q12-p22*q11+q11*q22-q12*q21),
             (p11*p22*q12-p11*q12*q22-p12*p21*q22-p12*p22*q11+p12*p22*q21+p12*q11*q22+p21*q12*q22-p22*q12*q21)/(p11*p22-p11*q22-p12*p21+p12*q21+p21*q12-p22*q11+q11*q22-q12*q21))

  int.pnt
} #end of the function
#'

#################################################################

#' @title The distance from a point to a line defined by two points
#'
#' @description Returns the distance from a point \code{p} to the line joining
#' points \code{a} and \code{b} in 2D space.
#'
#' @param p A 2D point, distance from \code{p} to the line passing through points \code{a} and \code{b} are to be
#' computed.
#' @param a,b 2D points that determine the straight line (i.e., through which the straight line passes).
#'
#' @return A \code{list} with two elements
#' \item{dis}{Distance from point \code{p} to the line passing through \code{a} and \code{b}}
#' \item{cl2p}{The closest point on the line passing through \code{a} and \code{b} to the point \code{p}}
#'
#' @seealso \code{\link{dist.point2plane}}, \code{\link{dist.point2set}} and \code{\link{Dist}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,2); B<-c(2,3); P<-c(3,1.5)
#'
#' dpl<-dist.point2line(P,A,B);
#' dpl
#' C<-dpl$cl2p
#' pts<-rbind(A,B,C,P)
#'
#' xr<-range(pts[,1])
#' xf<-(xr[2]-xr[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#' lnAB<-Line(A,B,x)
#' y<-lnAB$y
#' int<-lnAB$intercept #intercept
#' sl<-lnAB$slope #slope
#'
#' xsq<-seq(min(A[1],B[1],P[1])-xf,max(A[1],B[1],P[1])+xf,l=5)  #try also l=10, 20, or 100
#' pline<-(-1/sl)*(xsq-P[1])+P[2] #line passing thru P and perpendicular to AB
#'
#' Xlim<-range(pts[,1],x)
#' Ylim<-range(pts[,2],y)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(rbind(P),asp=1,pch=1,xlab="x",ylab="y",
#' main="Illustration of the distance from P \n to the Line Crossing Points A and B",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(rbind(A,B),pch=1)
#' lines(x,y,lty=1,xlim=Xlim,ylim=Ylim)
#' int<-my_round(int,2); sl<-my_round(sl,2)
#' text(rbind((A+B)/2+xd*c(-.01,-.01)),ifelse(sl==0,paste("y=",int),
#' ifelse(sl==1,paste("y=x+",int),
#' ifelse(int==0,paste("y=",sl,"x"),paste("y=",sl,"x+",int)))))
#' text(rbind(A+xd*c(0,-.01),B+xd*c(.0,-.01),P+xd*c(.01,-.01)),c("A","B","P"))
#' lines(xsq,pline,lty=2)
#' segments(P[1],P[2], C[1], C[2], lty=1,col=2,lwd=2)
#' text(rbind(C+xd*c(-.01,-.01)),"C")
#' text(rbind((P+C)/2),col=2,paste("d=",my_round(dpl$dis,2)))
#' }
#'
#' @export dist.point2line
dist.point2line <- function(p,a,b)
{
  if (!is.point(p) || !is.point(a) || !is.point(b))
  {stop('all arguments must be numeric points of dimension 2')}

  if (all(a==b))
  {stop('line is not well defined')}

  a1<-a[1]; a2<-a[2];
  b1<-b[1]; b2<-b[2];
  p1<-p[1]; p2<-p[2];

  if (a1==b1)
  {
    dp<-abs(a1-p1)
    C<-c(a1,p1)  #the point of intersection of perpendicular line from P to line joining A and B
  } else
  {
    sl<-slope(a,b);
    dp<-abs(p2-Line(a,b,p1)$y )/sqrt(sl^2+1)}
  c1<-(a1^2*p1-a1*a2*b2+a1*a2*p2-2*a1*b1*p1+a1*b2^2-a1*b2*p2+a2^2*b1-a2*b1*b2-a2*b1*p2+b1^2*p1+b1*b2*p2)/(a1^2-2*a1*b1+a2^2-2*a2*b2+b1^2+b2^2)
  c2<-(a1^2*b2-a1*a2*b1+a1*a2*p1-a1*b1*b2-a1*b2*p1+a2^2*p2+a2*b1^2-a2*b1*p1-2*a2*b2*p2+b1*b2*p1+b2^2*p2)/(a1^2-2*a1*b1+a2^2-2*a2*b2+b1^2+b2^2)
  C<-c(c1,c2)

  list(distance=dp, #dis is tge distance
       cl2p=C #cl2p closest point on the line to point p
  )


} #end of the functionP
#'

#################################################################

#' @title The distance between two vectors, matrices, or data frames
#'
#' @description Returns the Euclidean distance between \code{x} and \code{y} which can be vectors
#' or matrices or data frames of any dimension (\code{x} and \code{y} should be of same dimension).
#'
#' This function is different from the \code{\link[stats]{dist}} function in the \code{stats} package of the standard \code{R} distribution.
#' \code{dist} requires its argument to be a data matrix and \code{\link[stats]{dist}} computes and returns the distance matrix computed
#' by using the specified distance measure to compute the distances between the rows of a data matrix
#' (\insertCite{S-Book:1998;textual}{pcds}),
#' while \code{Dist} needs two arguments to find the distances between. For two data matrices \eqn{A} and \eqn{B},
#' \code{dist(rbind(as.vector(A),as.vector(B)))} and \code{Dist(A,B)} yield the same result.
#'
#' @param x,y Vectors, matrices or data frames (both should be of the same type).
#'
#' @return Euclidean distance between \code{x} and \code{y}
#'
#' @seealso \code{\link[stats]{dist}} from the base package \code{stats}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Dist(B,C);
#' dist(rbind(B,C))
#'
#' x<-runif(10)
#' y<-runif(10)
#' Dist(x,y)
#'
#' xm<-matrix(x,ncol=2)
#' ym<-matrix(y,ncol=2)
#' Dist(xm,ym)
#' dist(rbind(as.vector(xm),as.vector(ym)))
#'
#' Dist(xm,xm)
#' }
#'
#' @export Dist
Dist <- function(x,y)
{
  x<-as.matrix(x)
  y<-as.matrix(y)
  dis<-sqrt(sum((x-y)^2))

  dis
} #end of the function
#'

#################################################################

#' @title Distance from a point to a set of finite cardinality
#'
#' @description Returns the Euclidean distance between a point \code{p} and set of points \code{Yp} and the
#' closest point in set \code{Yp} to \code{p}. Distance between a point and a set is by definition the distance
#' from the point to the closest point in the set. \code{p} should be of finite dimension and \code{Yp} should
#' be of finite cardinality and \code{p} and elements of \code{Yp} must have the same dimension.
#'
#' @param p A \code{vector} (i.e., a point in \eqn{R^d}).
#' @param Yp A set of \eqn{d}-dimensional points.
#'
#' @return A \code{list} with the elements
#' \item{distance}{Distance from point \code{p} to set \code{Yp}}
#' \item{ind.cl.point}{Index of the closest point in set \code{Yp} to the point \code{p}}
#' \item{closest.point}{The closest point in set \code{Yp} to the point \code{p}}
#'
#' @seealso \code{\link{dist.point2line}} and \code{\link{dist.point2plane}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' dist.point2set(c(1,2),Te)
#'
#' X2<-cbind(runif(10),runif(10))
#' dist.point2set(c(1,2),X2)
#'
#' x<-runif(1)
#' y<-as.matrix(runif(10))
#' dist.point2set(x,y)  #this works, because x is a 1D point, and y is treated as a set of 10 1D points
#' #but will give an error message if y<-runif(10) is used above
#' }
#'
#' @export dist.point2set
dist.point2set <- function(p,Yp)
{
  if (!is.point(p,length(p)))
  {stop('p must be a numeric vector')}

  dim.p<-length(p)

  if (is.vector(Yp))
  {dim.Yp<-length(Yp);
  if (dim.Yp != dim.p )
  {stop('Both arguments must be of the same dimension')
  } else
  {dis<-Dist(p,Yp)
  ind.mdis<-1
  pcl<-Yp}
  } else
  {
    if (!is.matrix(Yp) && !is.data.frame(Yp))
    {stop('Yp must be a matrix or a data frame, each row representing a point')}
    Yp<-as.matrix(Yp)
    dim.Yp<-ncol(Yp);
    if (dim.Yp != dim.p)
    {stop('p and Yp must be of the same dimension')}
    dst<-vector()
    pcl<-NULL
    nY<-nrow(Yp)
    for (i in 1:nY)
    {
      dst<-c(dst,Dist(p,Yp[i,]) )
    }

    dis<-min(dst)
    ind.mdis<-which(dst==dis)
    pcl<-Yp[ind.mdis,]
  }

  list(
    distance=dis,
    ind.cl.point=ind.mdis,
    closest.point=pcl
  )
} #end of the function
#'

#################################################################

#' @title The function to round to \code{n} decimals or to the nearest non-zero decimal place, if rounding to \code{n} yields zero.
#'
#' @description
#' Returns the rounded number to the \code{n} decimals or to the nearest non-zero decimal place,
#' if rounding to \code{n} decimals yields zero. The default for \code{n} is 2.
#' Function is borrowed from this
#' \href{https://stackoverflow.com/}{link}.
#'
#' @param x A real number on which rounding will be applied
#' @param n Number of decimal places to round the number. If rounding yields 0, it rounds to the nearest non-zero decimal place.
#'
#' @return The rounded number to \code{n} decimals, and if this is zero to the nearest nonzero decimal
#'
#' @seealso \code{\link[base]{round}} and \code{\link[base]{signif}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' a=17; b=35
#'
#' my_round(a/b)
#' my_round(b/a)
#'
#' my_round(12.5624)
#' my_round(-0.64010)
#' my_round (0.000026)
#' my_round(-0.00072451)
#'
#' @export
my_round = function(x, n=2) {
  max(abs(round(x, n)), abs(signif(x, 1))) * sign(x)
} #end of the function
#'

#################################################################

#' @title The line at a point \code{p} parallel to the line segment joining two distinct 2D points \code{a} and \code{b}
#'
#' @description
#' An object of class \code{"Lines"}.
#' Returns the \code{equation, slope, intercept}, and \eqn{y}-coordinates of the line crossing the
#' point \code{p} and parallel to the line passing through the points \code{a} and \code{b} with
#' \eqn{x}-coordinates are provided in \code{vector} \code{x}.
#'
#' @param p A 2D point at which the parallel line to line segment joining \code{a} and \code{b} crosses.
#' @param a,b 2D points that determine the line segment (the line will be parallel to this line segment).
#' @param x A scalar or a \code{vector} of scalars representing the \eqn{x}-coordinates of the line parallel to
#' \code{ab} and crossing \code{p}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{Description of the line passing through point \code{p} and parallel to line segment joining
#' \code{a} and \code{b}}
#' \item{mtitle}{The \code{"main"} title for the plot of the line passing through point \code{p} and parallel to
#' line segment joining \code{a} and \code{b}}
#' \item{points}{The input points \code{p}, \code{a}, and \code{b} (stacked row-wise,
#' i.e., point \code{p} is in row 1, point \code{a} is in row 2 and point \code{b} is in row 3).
#' Line parallel to \code{ab} crosses p.}
#' \item{x}{The input vector. It can be a scalar or a \code{vector} of scalars,
#' which constitute the \eqn{x}-coordinates of the point(s) of interest on the line passing through point
#' \code{p} and parallel to line segment joining \code{a} and \code{b}.}
#' \item{y}{The output scalar or \code{vector} which constitutes the \eqn{y}-coordinates of the point(s) of interest
#' on the line passing through point \code{p} and parallel to line segment joining \code{a} and \code{b}.
#' If \code{x} is a scalar, then \code{y} will be a scalar and if \code{x} is a \code{vector} of scalars,
#' then \code{y} will be a \code{vector} of scalars.}
#' \item{slope}{Slope of the line, \code{Inf} is allowed, passing through point \code{p} and parallel to
#' line segment joining \code{a} and \code{b}}
#' \item{intercept}{Intercept of the line passing through point \code{p} and parallel to line segment
#' joining \code{a} and \code{b}}
#' \item{equation}{Equation of the line passing through point \code{p} and parallel to line segment joining
#' \code{a} and \code{b}}
#'
#' @seealso \code{\link{slope}}, \code{\link{Line}}, and \code{\link{perpline}}, \code{\link[stats]{line}}
#' in the generic \code{stats} package, \code{\link{paraline3D}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1.1,1.2); B<-c(2.3,3.4); p<-c(.51,2.5)
#'
#' paraline(p,A,B,.45)
#'
#' pts<-rbind(A,B,p)
#' xr<-range(pts[,1])
#' xf<-(xr[2]-xr[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#'
#' plnAB<-paraline(p,A,B,x)
#' plnAB
#' summary(plnAB)
#' plot(plnAB)
#'
#' y<-plnAB$y
#' Xlim<-range(x,pts[,1])
#' if (!is.na(y[1])) {Ylim<-range(y,pts[,2])} else {Ylim<-range(pts[,2])}
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' pf<-c(xd,-yd)*.025
#'
#' plot(A,pch=".",xlab="",ylab="",main="Line Crossing P and Parallel to AB",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(pts)
#' txt.str<-c("A","B","p")
#' text(pts+rbind(pf,pf,pf),txt.str)
#'
#' segments(A[1],A[2],B[1],B[2],lty=2)
#' if (!is.na(y[1])) {lines(x,y,type="l",lty=1,xlim=Xlim,ylim=Ylim)} else {abline(v=p[1])}
#' tx<-(A[1]+B[1])/2;
#' if (!is.na(y[1])) {ty<-paraline(p,A,B,tx)$y} else {ty=p[2]}
#' text(tx,ty,"line parallel to AB\n and crossing p")
#' }
#'
#' @export paraline
paraline <- function(p,a,b,x)
{
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))
  pname <-deparse(substitute(p))

  if (!is.point(p) || !is.point(a) || !is.point(b))
  {stop('p, a and b must all be numeric points of dimension 2')}

  if (!is.vector(x))
  {stop('x must be a vector')}

  if (all(a==b))
  {stop('a and b are same, hence lines are not well defined')}

  sl<-slope(a,b)
  if (abs(sl)==Inf)
  {
    int<-pln<-NA
  } else
  {
    int<-p[2]+sl*(0-p[1]);
    pln<-p[2]+sl*(x-p[1])
  }

  names(sl)<-"slope"
  names(int)<-"intercept"

  pts<-rbind(a,b,p)
  row.names(pts)<-c(aname,bname,pname)

  line.desc<-paste("Line Crossing Point ",pname, " and Parallel to Line Segment [",aname,bname,"]",sep="")
  main.text<-paste("Line Crossing Point ",pname, " and \n Parallel to Line Segment [",aname,bname,"]",sep="")

  if (abs(sl)==Inf)
  { vert<-p[1]
  eqn<-reqn<-paste("x = ",vert)
  } else
  {
    vert<-NA
    eqn<-ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),
                                            ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,
                                                                                                ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))

    rsl<-my_round(sl,2)
    rint<-my_round(int,2)
    reqn<-ifelse(rsl==0 & rint==0,"y=0",ifelse(rsl!=0 & rint==0,paste("y=",rsl,"x",sep=""),
                                               ifelse(rsl==0 & rint!=0,paste("y=",rint,sep=""),ifelse(rsl==1,
                                                                                                      ifelse(sign(rint)<0,paste("y=x",rint,sep=""),paste("y=x+",rint,sep="")),
                                                                                                      ifelse(sign(rint)<0,paste("y=",rsl,"x",rint,sep=""),paste("y=",rsl,"x+",rint,sep="")))) ))
  }

  res<-list(
    desc=line.desc, mtitle=main.text,
    points=pts, vert=vert,
    x=x,
    y=pln,
    slope=sl,
    intercept=int,
    equation=eqn,
    eqnlabel=reqn
  )
  class(res)<-"Lines"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Check whether a point is inside a circle
#'
#' @description Checks if the point \code{p} lies in the circle with center \code{cent} and radius \code{rad},
#' denoted as \code{C(cent,rad)}. So, it returns 1 or \code{TRUE} if \code{p} is inside the circle, and 0 otherwise.
#'
#' \code{boundary} is a logical argument (default=\code{FALSE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if the point, \code{p}, lies in the closure of the circle (i.e., interior and
#' boundary combined) else it checks if \code{p} lies in the interior of the circle.
#'
#' @param p A 2D point to be checked whether it is inside the circle or not.
#' @param cent A 2D point in Cartesian coordinates which serves as the center of the circle.
#' @param rad A positive real number which serves as the radius of the circle.
#' @param boundary A logical parameter (default=\code{TRUE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if the point, \code{p}, lies in the closure of the circle (i.e., interior and
#' boundary combined); else, it checks if \code{p} lies in the interior of the circle.
#'
#' @return Indicator for the point \code{p} being inside the circle or not, i.e., returns 1 or \code{TRUE}
#' if \code{p} is inside the circle, and 0 otherwise.
#'
#' @seealso \code{\link{in.triangle}}, \code{\link{in.tetrahedron}},  and
#' \code{\link[interp]{on.convex.hull}} from the \code{interp} package for documentation for \code{in.convex.hull}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' cent<-c(1,1); rad<-1; p<-c(1.4,1.2)
#' #try also cent<-runif(2); rad<-runif(1); p<-runif(2);
#'
#' in.circle(p,cent,rad)
#'
#' p<-c(.4,-.2)
#' in.circle(p,cent,rad)
#'
#' p<-c(1,0)
#' in.circle(p,cent,rad)
#' in.circle(p,cent,rad,boundary=FALSE)
#' }
#'
#' @export in.circle
in.circle <- function(p,cent,rad,boundary=TRUE)
{
  if (!is.point(p) || !is.point(cent) )
  {stop('p and cent must be numeric 2D points')}

  if (!is.point(rad,1) || rad<0)
  {stop('rad must be a scalar greater than or equal to 0')}

  if (boundary==TRUE)
  {ind.circ<-ifelse(Dist(p,cent)<=rad, TRUE, FALSE)
  }  else
  {ind.circ<-ifelse(Dist(p,cent)<rad, TRUE, FALSE)
  }
  ind.circ
} #end of the function
#'

#################################################################

#' @title The line passing through a point and perpendicular to the line segment joining two points
#'
#' @description
#' An object of class \code{"Lines"}.
#' Returns the \code{equation, slope, intercept}, and \eqn{y}-coordinates of the line crossing
#' the point \code{p} and perpendicular to the line passing through the points \code{a} and \code{b}
#' with \eqn{x}-coordinates are provided in \code{vector} \code{x}.
#'
#' @param p A 2D point at which the perpendicular line to line segment joining \code{a} and \code{b} crosses.
#' @param a,b 2D points that determine the line segment (the line will be perpendicular to this line segment).
#' @param x A scalar or a \code{vector} of scalars representing the \eqn{x}-coordinates of the line perpendicular to
#' line joining \code{a} and \code{b} and crossing \code{p}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{Description of the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}}
#' \item{mtitle}{The \code{"main"} title for the plot of the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}}
#' \item{points}{The input points \code{a} and \code{b} (stacked row-wise,
#' i.e., row 1 is point \code{a} and row 2 is point \code{b}). Line passing through point \code{p} is perpendicular
#' to line joining \code{a} and \code{b}}
#' \item{x}{The input vector, can be a scalar or a \code{vector} of scalars,
#' which constitute the \eqn{x}-coordinates of the point(s) of interest on the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}}
#' \item{y}{The output \code{vector} which constitutes the \eqn{y}-coordinates of the point(s) of interest on the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}. If \code{x} is a scalar, then \code{y} will be a scalar and if \code{x} is a \code{vector} of scalars,
#' then \code{y} will be a \code{vector} of scalars.}
#' \item{slope}{Slope of the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}}
#' \item{intercept}{Intercept of the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}}
#' \item{equation}{Equation of the line passing through point \code{p} and perpendicular to line joining
#' \code{a} and \code{b}}
#'
#' @seealso \code{\link{slope}}, \code{\link{Line}}, and \code{\link{paraline}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1.1,1.2); B<-c(2.3,3.4); p<-c(.51,2.5)
#'
#' perpline(p,A,B,.45)
#'
#' pts<-rbind(A,B,p)
#' xr<-range(pts[,1])
#' xf<-(xr[2]-xr[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#'
#' plnAB<-perpline(p,A,B,x)
#' plnAB
#' summary(plnAB)
#' plot(plnAB,asp=1)
#'
#' y<-plnAB$y
#' Xlim<-range(x,pts[,1])
#' if (!is.na(y[1])) {Ylim<-range(y,pts[,2])} else {Ylim<-range(pts[,2])}
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' pf<-c(xd,-yd)*.025
#'
#' plot(A,asp=1,pch=".",xlab="",ylab="",main="Line Crossing p and Perpendicular to AB",
#' xlim=Xlim+xd*c(-.5,.5),ylim=Ylim+yd*c(-.05,.05))
#' points(pts)
#' txt.str<-c("A","B","p")
#' text(pts+rbind(pf,pf,pf),txt.str)
#'
#' segments(A[1],A[2],B[1],B[2],lty=2)
#' if (!is.na(y[1])) {lines(x,y,type="l",lty=1,xlim=Xlim,ylim=Ylim)} else {abline(v=p[1])}
#' tx<-p[1]+abs(xf-p[1])/2;
#' if (!is.na(y[1])) {ty<-perpline(p,A,B,tx)$y} else {ty=p[2]}
#' text(tx,ty,"line perpendicular to AB\n and crossing p")
#' }
#'
#' @export perpline
perpline <- function(p,a,b,x)
{
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))
  pname <-deparse(substitute(p))

  if (!is.point(p) || !is.point(a) || !is.point(b))
  {stop('p, a and b must all be numeric points of dimension 2')}

  if (all(a==b))
  {stop('a and b are same, hence lines are not well defined')}

  if (!is.vector(x))
  {stop('x must be a vector')}

  sl<--1/slope(a,b);
  int<-p[2]+sl *(0-p[1]);
  perpln<-p[2]+sl *(x-p[1])

  if (abs(sl)==Inf)
  {
    int<-perpln<-NA
  } else
  {
    int<-p[2]+sl *(0-p[1]);
    perpln<-p[2]+sl *(x-p[1])
  }

  names(sl)<-"slope"
  names(int)<-"intercept"

  line.desc<-paste("Line Crossing Point ",pname, " Perpendicular to Line Segment [",aname,bname,"]",sep="")
  main.txt<-paste("Line Crossing Point ",pname, " and \n Perpendicular to Line Segment [",aname,bname,"]",sep="")

  pts<-rbind(a,b,p)
  row.names(pts)<-c(aname,bname,pname)

  if (abs(sl)==Inf)
  { vert<-p[1]
  eqn<-reqn<-paste("x =",vert)
  } else
  {
    vert<-NA
    eqn<-ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),
                                            ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,
                                                                                                ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
    rsl<-my_round(sl,2)
    rint<-my_round(int,2)
    reqn<-ifelse(rsl==0 & rint==0,"y=0",ifelse(rsl!=0 & rint==0,paste("y=",rsl,"x",sep=""),
                                               ifelse(rsl==0 & rint!=0,paste("y=",rint,sep=""),ifelse(rsl==1,
                                                                                                      ifelse(sign(rint)<0,paste("y=x",rint,sep=""),paste("y=x+",rint,sep="")),
                                                                                                      ifelse(sign(rint)<0,paste("y=",rsl,"x",rint,sep=""),paste("y=",rsl,"x+",rint,sep="")))) ))
  }

  res<-list(
    desc=line.desc, mtitle=main.txt,
    points=pts, vert=vert,
    x=x,
    y=perpln,
    slope=sl,
    intercept=int,
    equation=eqn,
    eqnlabel=reqn
  )

  class(res)<-"Lines"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The radii of points from one class with respect to points from the other class
#'
#' @description Returns the radii of the balls centered at \code{x} points where radius of an \code{x} point equals to the minimum distance
#' to \code{y} points (i.e., distance to the closest \code{y} point). That is, for each \code{x} point \eqn{radius= \min_{y \in Y}(d(x,y))}.
#' \code{x} and \code{y} points must be of the same dimension.
#'
#' @param x A set of \eqn{d}-dimensional points for which the radii are computed. Radius of an \code{x} point equals to the
#' distance to the closest \code{y} point.
#' @param y A set of \eqn{d}-dimensional points representing the reference points for the balls. That is, radius
#' of an \code{x} point is defined as the minimum distance to the \code{y} points.
#'
#' @return A \code{list} with three elements
#' \item{rad}{A \code{vector} whose entries are the radius values for the \code{x} points. Radius of an \code{x} point equals to
#' the distance to the closest \code{y} point}
#' \item{index.of.clYp}{A \code{vector} of indices of the closest \code{y} points to the \code{x} points. The \eqn{i}-th entry in this
#' \code{vector} is the index of the closest \code{y} point to \eqn{i}-th \code{x} point.}
#' \item{closest.Yp}{A \code{vector} of the closest \code{y} points to the \code{x} points. The \eqn{i}-th entry in this
#' \code{vector} or \eqn{i}-th row in the matrix is the closest \code{y} point to \eqn{i}-th \code{x} point.}
#'
#' @seealso \code{\link{radius}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-10
#' ny<-5
#' X<-cbind(runif(nx),runif(nx))
#' Y<-cbind(runif(ny),runif(ny))
#' Rad<-radii(X,Y)
#' Rad
#' rd<-Rad$rad
#'
#' Xlim<-range(X[,1]-rd,X[,1]+rd,Y[,1])
#' Ylim<-range(X[,2]-rd,X[,2]+rd,Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(rbind(Y),asp=1,pch=16,col=2,xlab="",ylab="",
#' main="Circles Centered at Class X Points with \n Radius Equal to the Distance to Closest Y Point",
#' axes=TRUE, xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(rbind(X))
#' interp::circles(X[,1],X[,2],Rad$rad,lty=1,lwd=1,col=4)
#'
#' #For 1D data
#' nx<-10
#' ny<-5
#' Xm<-as.matrix(X)
#' Ym<-as.matrix(Y)
#' radii(Xm,Ym)  #this works as Xm and Ym are treated as 1D data sets
#' #but will give error if radii(X,Y) is used
#' #as X and Y are treated as vectors (i.e., points)
#'
#' #For 3D data
#' nx<-10
#' ny<-5
#' X<-cbind(runif(nx),runif(nx),runif(nx))
#' Y<-cbind(runif(ny),runif(ny),runif(ny))
#' radii(X,Y)
#' }
#'
#' @export radii
radii <- function(x,y)
{
  if (!is.numeric(as.matrix(x)) || !is.numeric(as.matrix(y)) )
  {stop('both arguments must be numeric')}

  ifelse(is.point(x,length(x)), x<-matrix(x,nrow=1),x<-as.matrix(x))
  ifelse(is.point(y,length(y)), y<-matrix(y,nrow=1),y<-as.matrix(y))

  dimx<-dimension(x)
  dimy<-dimension(y)

  if (dimx != dimy)
  {stop('the arguments are not of the same dimension')}

  nx<-nrow(x)
  ny<-nrow(y)

  rad<-rep(0,nx)
  cl.ind<-vector()
  for (i in 1:nx)
  { dis<-rep(0,ny)
  for (j in 1:ny)
  {
    dis[j]<-Dist(x[i,],y[j,])
  }
  rad[i]<-min(dis)
  cl.ind<-c(cl.ind,which(dis==rad[i]))
  }
  list(
    radiuses=rad,
    indices.of.closest.points=cl.ind,
    closest.points=y[cl.ind,]
  )
} #end of the function
#'

#################################################################

#' @title The radius of a point from one class with respect to points from the other class
#'
#' @description Returns the radius for the ball centered at point \code{p} with radius=min distance to \code{Y} points.
#' That is, for the point \code{p} \eqn{radius= \min_{y \in Y}d(p,y)} (i.e., distance from \code{p} to the closest \code{Y} point).
#' The point \code{p} and \code{Y} points must be of same dimension.
#'
#' @param p A \eqn{d}-dimensional point for which radius is computed. Radius of \code{p} equals to the
#' distance to the closest \code{Y} point to \code{p}.
#' @param Y A set of \eqn{d}-dimensional points representing the reference points for the balls. That is, radius
#' of the point \code{p} is defined as the minimum distance to the \code{Y} points.
#'
#' @return A \code{list} with three elements
#' \item{rad}{Radius value for the point, \code{p} defined as \eqn{\min_{y in Y} d(p,y)}}
#' \item{index.of.clYpnt}{Index of the closest \code{Y} points to the point \code{p}}
#' \item{closest.Ypnt}{The closest \code{Y} point to the point \code{p}}
#'
#' @seealso \code{\link{radii}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#'
#' ny<-10
#' Y<-cbind(runif(ny),runif(ny))
#' radius(A,Y)
#'
#' nx<-10
#' X<-cbind(runif(nx),runif(nx))
#' rad<-rep(0,nx)
#' for (i in 1:nx)
#' rad[i]<-radius(X[i,],Y)$rad
#'
#' Xlim<-range(X[,1]-rad,X[,1]+rad,Y[,1])
#' Ylim<-range(X[,2]-rad,X[,2]+rad,Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(rbind(Y),asp=1,pch=16,col=2,xlab="",ylab="",
#' main="Circles Centered at Class X Points with \n Radius Equal to the Distance to Closest Y Point",
#' axes=TRUE, xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(rbind(X))
#' interp::circles(X[,1],X[,2],rad,lty=1,lwd=1,col=4)
#'
#' #For 1D data
#' ny<-5
#' Y<-runif(ny)
#' Ym = as.matrix(Y)
#' radius(1,Ym)  #this works as Y is treated as 1D data sets
#' #but will give error if radius(1,Y) is used
#' #as Y is treated as a vector (i.e., points)
#'
#' #For 3D data
#' ny<-5
#' X<-runif(3)
#' Y<-cbind(runif(ny),runif(ny),runif(ny))
#' radius(X,Y)
#' }
#'
#' @export radius
radius <- function(p,Y)
{
  if (!is.point(p,length(p)))
  {stop('p must be a numeric point')}

  if (!is.numeric(as.matrix(Y)) )
  {stop('Y must be numeric')}

  ifelse(is.point(Y,length(Y)), Y<-matrix(Y,nrow=1),Y<-as.matrix(Y))

  ny<-nrow(Y)

  dis<-rep(0,ny)
  for (j in 1:ny)
    dis[j]<-Dist(p,Y[j,]);

  min.dis<-min(dis)
  cl.ind<-which(dis==min.dis)
  cl.Ypnt<-Y[cl.ind,]
  list(
    rad=min.dis,
    index.of.clYpnt=cl.ind,
    closest.Ypnt=cl.Ypnt
  )
} #end of the function
#'

#################################################################

#' @title The area of a polygon in \eqn{R^2}
#'
#' @description Returns the area of the polygon, \code{h}, in the real plane \eqn{R^{2}}; the vertices of the polygon \code{h}
#' must be provided in clockwise or counter-clockwise order, otherwise the function does not yield
#' the area of the polygon. Also, the polygon could be convex or non-convex.
#' See (\insertCite{weisstein-area-poly;textual}{pcds}).
#'
#' @param h A \code{vector} of \eqn{n} 2D points, stacked row-wise, each row representing a vertex of the polygon,
#' where \eqn{n} is the number of vertices of the polygon.
#'
#' @return area of the polygon \code{h}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(0.5,.8);
#' Tr<-rbind(A,B,C);
#' area.polygon(Tr)
#'
#' A<-c(0,0); B<-c(1,0); C<-c(.7,.6); D<-c(0.3,.8);
#' h1<-rbind(A,B,C,D);  #try also h1<-rbind(A,B,D,C) or h1<-rbind(A,C,B,D) or h1<-rbind(A,D,C,B);
#' area.polygon(h1)
#'
#' Xlim<-range(h1[,1])
#' Ylim<-range(h1[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(h1,xlab="",ylab="",main="A Convex Polygon with Four Vertices",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(h1)
#' xc<-rbind(A,B,C,D)[,1]+c(-.03,.03,.02,-.01)
#' yc<-rbind(A,B,C,D)[,2]+c(.02,.02,.02,.03)
#' txt.str<-c("A","B","C","D")
#' text(xc,yc,txt.str)
#'
#' #when the triangle is degenerate, it gives zero area
#' B<-A+2*(C-A);
#' T2<-rbind(A,B,C)
#' area.polygon(T2)
#' }
#'
#' @export area.polygon
area.polygon <- function(h)
{
  h<-as.matrix(h)

  if (!is.numeric(h) || ncol(h)!=2)
  {stop('the argument must be numeric and of dimension np2')}

  n<-nrow(h);
  area<-0;
  for (i in 1:(n-1) )
  { a1<-h[i,1]; b1<-h[i,2]; a2<-h[i+1,1]; b2<-h[i+1,2];
  area<-area + 1/2*(a1*b2-a2*b1);
  }
  a1<-h[1,1]; b1<-h[1,2];

  abs(as.numeric(area + 1/2*(a2*b1-a1*b2)));
} #end of the function
#'

#################################################################

#' @title The plane at a point and parallel to the plane spanned by three distinct 3D points \code{a}, \code{b}, and \code{c}
#'
#' @description
#' An object of class \code{"Planes"}.
#' Returns the equation and \eqn{z}-coordinates of the plane passing through point \code{p} and parallel to the plane spanned
#' by three distinct 3D points \code{a}, \code{b}, and \code{c} with \eqn{x}- and \eqn{y}-coordinates are provided in vectors \code{x} and \code{y},
#' respectively.
#'
#' @param p A 3D point which the plane parallel to the plane spanned by
#' three distinct 3D points \code{a}, \code{b}, and \code{c} crosses.
#' @param a,b,c 3D points that determine the plane to which the plane crossing point \code{p} is parallel to.
#' @param x,y Scalars or \code{vectors} of scalars representing the \eqn{x}- and \eqn{y}-coordinates of the plane parallel to
#' the plane spanned by points \code{a}, \code{b}, and \code{c} and passing through point \code{p}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{Description of the plane passing through point \code{p} and parallel to plane spanned by points
#' \code{a}, \code{b} and \code{c}}
#' \item{points}{The input points \code{a}, \code{b}, \code{c}, and \code{p}.
#' Plane is parallel to the plane spanned by \code{a}, \code{b}, and \code{c}
#' and passes through point \code{p} (stacked row-wise, i.e., row 1 is point \code{a}, row 2 is point \code{b},
#'  row 3 is point \code{c}, and row 4 is point \code{p}).}
#' \item{x,y}{The input vectors which constitutes the \eqn{x}- and \eqn{y}-coordinates of the point(s) of interest on the
#' plane. \code{x} and \code{y} can be scalars or vectors of scalars.}
#' \item{z}{The output \code{vector} which constitutes the \eqn{z}-coordinates of the point(s) of interest on the plane.
#' If \code{x} and \code{y} are scalars, \code{z} will be a scalar and
#' if \code{x} and \code{y} are vectors of scalars, then \code{z} needs to be a \code{matrix} of scalars,
#' containing the \eqn{z}-coordinate for each pair of \code{x} and \code{y} values.}
#' \item{coeff}{Coefficients of the plane (in the \eqn{z = A x+B y+C} form).}
#' \item{equation}{Equation of the plane in long form}
#' \item{equation2}{Equation of the plane in short form, to be inserted on the plot}
#'
#' @seealso \code{\link{Plane}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' Q<-c(1,10,3); R<-c(1,1,3); S<-c(3,9,12); P<-c(1,1,0)
#'
#' pts<-rbind(Q,R,S,P)
#' paraplane(P,Q,R,S,.1,.2)
#'
#' xr<-range(pts[,1]); yr<-range(pts[,2])
#' xf<-(xr[2]-xr[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' yf<-(yr[2]-yr[1])*.25 #how far to go at the lower and upper ends in the y-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#' y<-seq(yr[1]-yf,yr[2]+yf,l=5)  #try also l=10, 20, or 100
#'
#' plP2QRS<-paraplane(P,Q,R,S,x,y)
#' plP2QRS
#' summary(plP2QRS)
#' plot(plP2QRS,theta = 225, phi = 30, expand = 0.7, facets = FALSE, scale = TRUE)
#'
#' paraplane(P,Q,R,Q+R,.1,.2)
#'
#' z.grid<-plP2QRS$z
#'
#' plQRS<-Plane(Q,R,S,x,y)
#' plQRS
#' pl.grid<-plQRS$z
#'
#' zr<-max(z.grid)-min(z.grid)
#' Pts<-rbind(Q,R,S,P)+rbind(c(0,0,zr*.1),c(0,0,zr*.1),c(0,0,zr*.1),c(0,0,zr*.1))
#' Mn.pts<-apply(Pts[1:3,],2,mean)
#'
#' plot3D::persp3D(z = pl.grid, x = x, y = y, theta =225, phi = 30, ticktype = "detailed",
#' main="Plane Crossing Points Q, R, S\n and Plane Passing P Parallel to it")
#' #plane spanned by points Q, R, S
#' plot3D::persp3D(z = z.grid, x = x, y = y,add=TRUE)
#' #plane parallel to the original plane and passing thru point \code{P}
#'
#' plot3D::persp3D(z = z.grid, x = x, y = y, theta =225, phi = 30, ticktype = "detailed",
#' main="Plane Crossing Point P \n and Parallel to the Plane Crossing Q, R, S")
#' #plane spanned by points Q, R, S
#' #add the defining points
#' plot3D::points3D(Pts[,1],Pts[,2],Pts[,3], add=TRUE)
#' plot3D::text3D(Pts[,1],Pts[,2],Pts[,3], c("Q","R","S","P"),add=TRUE)
#' plot3D::text3D(Mn.pts[1],Mn.pts[2],Mn.pts[3],plP2QRS$equation,add=TRUE)
#' plot3D::polygon3D(Pts[1:3,1],Pts[1:3,2],Pts[1:3,3], add=TRUE)
#' }
#'
#' @export paraplane
paraplane <- function(p,a,b,c,x,y)
{
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))
  cname <-deparse(substitute(c))
  pname <-deparse(substitute(p))

  if (!is.point(p,3) || !is.point(a,3) || !is.point(b,3) || !is.point(c,3))
  {stop('p, a, b, and c must be numeric points of dimension 3')}

  dab<-Dist(a,b); dac<-Dist(a,c); dbc<-Dist(b,c);
  sd<-sort(c(dab,dac,dbc))

  if (isTRUE(all.equal(sd[3] , sd[1]+sd[2])))
  {stop('the three defining points, a, b, and c are collinear; plane is not well-defined')}

  if (!is.point(x,length(x)) || !is.point(y,length(y)))
  {stop('x and y must be numeric vectors')}

  p1<-p[1]; p2<-p[2]; p3<-p[3];
  a1<-a[1]; a2<-a[2]; a3<-a[3];
  b1<-b[1]; b2<-b[2]; b3<-b[3];
  c1<-c[1]; c2<-c[2]; c3<-c[3];

  A <- (-a2*b3+a2*c3+a3*b2-a3*c2-b2*c3+b3*c2)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1)
  B <- (a1*b3-a1*c3-a3*b1+a3*c1+b1*c3-b3*c1)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1)
  C <- (a1*b2*p3-a1*b3*p2-a1*c2*p3+a1*c3*p2-a2*b1*p3+a2*b3*p1+a2*c1*p3-a2*c3*p1+a3*b1*p2-a3*b2*p1-a3*c1*p2+a3*c2*p1+b1*c2*p3-b1*c3*p2-b2*c1*p3+b2*c3*p1+b3*c1*p2-b3*c2*p1)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1)

  z.grid <- outer(x, y, function(a, b) A*a+B*b+C)

  plane.desc<-paste("Plane Passing through Point ", pname," Parallel to the Plane \n Passing through Points ", aname,", ",bname, " and ", cname,sep="")
  plane.title<-paste("Plane at Point ", pname," Parallel to the Plane \n Crossing Points ", aname,", ",bname, " and ", cname,sep="")

  pts<-rbind(a,b,c,p)
  row.names(pts)<-c(aname,bname,cname,pname)

  condA= c(A==0,A==1,!(A==0|A==1))
  condB= c(B==0,B==1,!(B==0|B==1))
  compx =switch(which(condA==TRUE),"",paste("x",sep=""),paste(A,"x",sep=""))
  compy =switch(which(condB==TRUE),"",paste("y",sep=""),paste(B,"y",sep=""))

  signB=ifelse(sign(B)<=0,""," + ")
  signC=ifelse(sign(C)<=0," "," + ")

  eqn=paste("z = ",compx,signB,compy,signC,C,sep="");

  compp2 =switch(which(condA==TRUE),"",paste("x",sep=""),paste(my_round(A,2),"x",sep=""))
  compy2 =switch(which(condB==TRUE),"",paste("y",sep=""),paste(my_round(B,2),"y",sep=""))
  eqn2=paste("z = ",compp2,signB,compy2,signC,my_round(C,2),sep="");

  coeffs<-c(A,B,C)
  names(coeffs)<-c("A","B","C")

  res<-list(
    desc=plane.desc,
    main.title=plane.title,
    points=pts,
    x=x,
    y=y,
    z=z.grid,
    coeff=coeffs,
    equation=eqn,
    equation2=eqn2
  )

  class(res)<-"Planes"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################
#Functions in R^3
#################################################################

#' @title The line crossing 3D point \code{p} in the direction of \code{vector} \code{v} (or if \code{v} is a point,
#' in direction of \eqn{v-r_0})
#'
#' @description
#' An object of class \code{"Lines3D"}.
#' Returns the equation, \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinates of the line crossing 3D point \eqn{r_0}
#' in the direction of \code{vector} \code{v} (of if \code{v} is a point, in the direction of \eqn{v-r_0})
#' with the parameter \code{t} being provided in \code{vector} \code{t}.
#'
#' @param p A 3D point through which the straight line passes.
#' @param v A 3D \code{vector} which determines the direction of the straight line (i.e., the straight line would be
#' parallel to this vector) if the \code{dir.vec=TRUE}, otherwise it is 3D point and \eqn{v-r_0} determines the direction of the
#' the straight line.
#' @param t A scalar or a \code{vector} of scalars representing the parameter of the coordinates of the line
#' (for the form: \eqn{x=p_0 + a t}, \eqn{y=y_0 + b t}, and \eqn{z=z_0 + c t} where \eqn{r_0=(p_0,y_0,z_0)}
#' and \eqn{v=(a,b,c)} if \code{dir.vec=TRUE}, else \eqn{v-r_0=(a,b,c)}).
#' @param dir.vec A logical argument about \code{v}, if \code{TRUE} \code{v} is treated as a vector,
#' else \code{v} is treated as a point and so the direction \code{vector} is taken to be \eqn{v-r_0}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A description of the line}
#' \item{mtitle}{The \code{"main"} title for the plot of the line}
#' \item{pts}{The input points that determine a line and/or a plane, \code{NULL} for this function.}
#' \item{pnames}{The names of the input points that determine a line and/or a plane, \code{NULL} for this function.}
#' \item{vecs}{The point \code{p} and the \code{vector} \code{v} (if \code{dir.vec=TRUE}) or the point \code{v} (if \code{dir.vec=FALSE}). The first row is \code{p}
#' and the second row is \code{v}.}
#' \item{vec.names}{The names of the point \code{p} and the \code{vector} \code{v} (if \code{dir.vec=TRUE}) or the point \code{v} (if \code{dir.vec=FALSE}).}
#' \item{x,y,z}{The \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinates of the point(s) of interest on the line.}
#' \item{tsq}{The scalar or the \code{vector} of the parameter in defining each coordinate of the line for the form:
#' \eqn{x=p_0 + a t}, \eqn{y=y_0 + b t}, and \eqn{z=z_0 + c t} where \eqn{r_0=(p_0,y_0,z_0)} and \eqn{v=(a,b,c)} if \code{dir.vec=TRUE}, else \eqn{v-r_0=(a,b,c)}.}
#' \item{equation}{Equation of the line passing through point \code{p} in the direction of the \code{vector} \code{v} (if \code{dir.vec=TRUE})
#' else in the direction of \eqn{v-r_0}. The line equation is in the form: \eqn{x=p_0 + a t}, \eqn{y=y_0 + b t}, and \eqn{z=z_0 + c t} where
#' \eqn{r_0=(p_0,y_0,z_0)} and \eqn{v=(a,b,c)} if \code{dir.vec=TRUE}, else \eqn{v-r_0=(a,b,c)}.}
#'
#' @seealso \code{\link{line}}, \code{\link{paraline3D}}, and \code{\link{Plane}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,10,3); B<-c(1,1,3);
#'
#' vecs<-rbind(A,B)
#'
#' Line3D(A,B,.1)
#' Line3D(A,B,.1,dir.vec=FALSE)
#'
#' tr<-range(vecs);
#' tf<-(tr[2]-tr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' tsq<-seq(-tf*10-tf,tf*10+tf,l=5)  #try also l=10, 20, or 100
#'
#' lnAB3D<-Line3D(A,B,tsq)  #try also lnAB3D<-Line3D(A,B,tsq,dir.vec=FALSE)
#' lnAB3D
#' summary(lnAB3D)
#' plot(lnAB3D)
#'
#' x<-lnAB3D$x
#' y<-lnAB3D$y
#' z<-lnAB3D$z
#'
#' zr<-range(z)
#' zf<-(zr[2]-zr[1])*.2
#' Bv<-B*tf*5
#'
#' Xlim<-range(x)
#' Ylim<-range(y)
#' Zlim<-range(z)
#'
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' Dr<-A+min(tsq)*B
#'
#' plot3D::lines3D(x, y, z, phi = 0, bty = "g",main="Line Crossing A \n in the Direction of OB",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.1,.1),
#'         pch = 20, cex = 2, ticktype = "detailed")
#' plot3D::arrows3D(Dr[1],Dr[2],Dr[3]+zf,Dr[1]+Bv[1],Dr[2]+Bv[2],Dr[3]+zf+Bv[3], add=TRUE)
#' plot3D::points3D(A[1],A[2],A[3],add=TRUE)
#' plot3D::arrows3D(A[1],A[2],A[3]-2*zf,A[1],A[2],A[3],lty=2, add=TRUE)
#' plot3D::text3D(A[1],A[2],A[3]-2*zf,labels="initial point",add=TRUE)
#' plot3D::text3D(A[1],A[2],A[3]+zf/2,labels=expression(r[0]),add=TRUE)
#' plot3D::arrows3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+3*zf+Bv[3]/2,
#' Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+zf+Bv[3]/2,lty=2, add=TRUE)
#' plot3D::text3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+3*zf+Bv[3]/2,
#' labels="direction vector",add=TRUE)
#' plot3D::text3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+zf+Bv[3]/2,labels="v",add=TRUE)
#' plot3D::text3D(0,0,0,labels="O",add=TRUE)
#' }
#'
#' @export Line3D
Line3D <- function(p,v,t,dir.vec=TRUE)
{
  pname <-deparse(substitute(p))
  vname <-ifelse(dir.vec,paste("O",deparse(substitute(v)),sep=""),
                 paste(pname,deparse(substitute(v)),sep=""))

  if (!is.point(p,3) || !is.point(v,3))
  {stop('p and v must be numeric vectors of dimension 3')}

  if (!is.point(t,length(t)))
  {stop('t must be a numeric vector')}

  if (isTRUE(all.equal(v,c(0,0,0))))
  {stop('The line is not well-defined, as the direction vector is the zero vector!')}

  x0<-p[1]; y0<-p[2]; z0<-p[3];

  if (!dir.vec) {v <- v-p}

  a<-v[1]; b<-v[2]; c<-v[3];

  x=x0+a*t
  y=y0+b*t
  z=z0+c*t

  vecs<-rbind(p,v)
  row.names(vecs) <- c("initial point","direction vector")

  ifelse(dir.vec,
  line.desc<-paste("the line passing through point", pname, "in the direction of", vname, " with O representing the origin (0,0,0)
                   (i.e., parallel to", vname, ")"),
  line.desc<-paste("the line passing through point", pname, "in the direction of", vname, "(i.e., parallel to", vname, ")"))
  main.txt<-paste("Line Crossing Point", pname, "\n in the Direction of", vname)

  conda= c(a==0,a==1,a==-1,!(a==0|abs(a)==1))
  condb= c(b==0,b==1,b==-1,!(b==0|abs(b)==1))
  condc= c(c==0,c==1,c==-1,!(c==0|abs(c)==1))
  cx =switch(which(conda==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(a,"t",sep=""))
  cy =switch(which(condb==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(b,"t",sep=""))
  cz =switch(which(condc==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(c,"t",sep=""))
  cx0 =ifelse(x0==0,"",x0)
  cy0 =ifelse(y0==0,"",y0)
  cz0 =ifelse(z0==0,"",z0)
  signa=ifelse(sign(a)<=0,""," + ")
  signb=ifelse(sign(b)<=0,""," + ")
  signc=ifelse(sign(c)<=0,""," + ")

  eqn<- rbind(paste("x = ",cx0,signa,cx,sep=""),
              paste("y = ",cy0,signb,cy,sep=""),
              paste("z = ",cz0,signc,cz,sep=""))

  res<-list(
    desc=line.desc,
    mtitle=main.txt,
    pts=NULL,
    pnames=NULL,
    vecs=vecs,
    vec.names=c(pname,vname),
    x=x,
    y=y,
    z=z,
    tsq=t,
    equation=eqn
  )
  class(res)<-"Lines3D"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The line crossing the 3D point \code{p} and parallel to line joining 3D points \code{a} and \code{b}
#'
#' @description
#' An object of class \code{"Lines3D"}.
#' Returns the equation, \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinates of the line crossing 3D point \code{p} and parallel to the line
#' joining 3D points \code{a} and \code{b} (i.e., the line is in the direction of \code{vector} \code{b}-\code{a})
#' with the parameter \code{t} being provided in \code{vector} \code{t}.
#'
#' @param p A 3D point through which the straight line passes.
#' @param a,b 3D points which determine the straight line to which the line passing through point \code{p} would be
#' parallel (i.e., \eqn{b-a} determines the direction of the straight line passing through \code{p}).
#' @param t A scalar or a \code{vector} of scalars representing the parameter of the coordinates of the line
#' (for the form: \eqn{x=p_0 + A t}, \eqn{y=y_0 + B t}, and \eqn{z=z_0 + C t} where \eqn{p=(p_0,y_0,z_0)}
#' and \eqn{b-a=(A,B,C)}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A description of the line}
#' \item{mtitle}{The \code{"main"} title for the plot of the line}
#' \item{points}{The input points that determine the line to which the line crossing point \code{p} would be parallel.}
#' \item{pnames}{The names of the input points that determine the line to which the line crossing point \code{p} would
#' be parallel.}
#' \item{vecs}{The points \code{p}, \code{a}, and \code{b} stacked row-wise in this order.}
#' \item{vec.names}{The names of the points \code{p}, \code{a}, and \code{b}.}
#' \item{x,y,z}{The \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinates of the point(s) of interest on the line parallel to the line
#' determined by points \code{a} and \code{b}.}
#' \item{tsq}{The scalar or the \code{vector} of the parameter in defining each coordinate of the line for the form:
#' \eqn{x=p_0 + A t}, \eqn{y=y_0 + B t}, and \eqn{z=z_0 + C t} where \eqn{p=(p_0,y_0,z_0)} and \eqn{b-a=(A,B,C)}.}
#' \item{equation}{Equation of the line passing through point \code{p} and parallel to the line joining points \code{a} and \code{b}
#' (i.e., in the direction of the \code{vector} \code{b}-\code{a}). The line equation is in the form: \eqn{x=p_0 + A t}, \eqn{y=y_0 + B t}, and
#' \eqn{z=z_0 + C t} where \eqn{p=(p_0,y_0,z_0)} and \eqn{b-a=(A,B,C)}.}
#'
#' @seealso \code{\link{Line3D}}, \code{\link{perpline2plane}}, and \code{\link{paraline}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P<-c(1,10,4); Q<-c(1,1,3); R<-c(3,9,12)
#'
#' vecs<-rbind(P,R-Q)
#' pts<-rbind(P,Q,R)
#' paraline3D(P,Q,R,.1)
#'
#' tr<-range(pts,vecs);
#' tf<-(tr[2]-tr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' tsq<-seq(-tf*10-tf,tf*10+tf,l=5)  #try also l=10, 20, or 100
#'
#' pln3D<-paraline3D(P,Q,R,tsq)
#' pln3D
#' summary(pln3D)
#' plot(pln3D)
#'
#' x<-pln3D$x
#' y<-pln3D$y
#' z<-pln3D$z
#'
#' zr<-range(z)
#' zf<-(zr[2]-zr[1])*.2
#' Qv<-(R-Q)*tf*5
#'
#' Xlim<-range(x,pts[,1])
#' Ylim<-range(y,pts[,2])
#' Zlim<-range(z,pts[,3])
#'
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' Dr<-P+min(tsq)*(R-Q)
#'
#' plot3D::lines3D(x, y, z, phi = 0, bty = "g", main="Line Crossing P \n in the direction of R-Q",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.1,.1)+c(-zf,zf),
#'         pch = 20, cex = 2, ticktype = "detailed")
#' plot3D::arrows3D(Dr[1],Dr[2],Dr[3]+zf,Dr[1]+Qv[1],Dr[2]+Qv[2],Dr[3]+zf+Qv[3], add=TRUE)
#' plot3D::points3D(pts[,1],pts[,2],pts[,3],add=TRUE)
#' plot3D::text3D(pts[,1],pts[,2],pts[,3],labels=c("P","Q","R"),add=TRUE)
#' plot3D::arrows3D(P[1],P[2],P[3]-2*zf,P[1],P[2],P[3],lty=2, add=TRUE)
#' plot3D::text3D(P[1],P[2],P[3]-2*zf,labels="initial point",add=TRUE)
#' plot3D::arrows3D(Dr[1]+Qv[1]/2,Dr[2]+Qv[2]/2,Dr[3]+3*zf+Qv[3]/2,Dr[1]+Qv[1]/2,
#' Dr[2]+Qv[2]/2,Dr[3]+zf+Qv[3]/2,lty=2, add=TRUE)
#' plot3D::text3D(Dr[1]+Qv[1]/2,Dr[2]+Qv[2]/2,Dr[3]+3*zf+Qv[3]/2,labels="direction vector",add=TRUE)
#' plot3D::text3D(Dr[1]+Qv[1]/2,Dr[2]+Qv[2]/2,Dr[3]+zf+Qv[3]/2,labels="R-Q",add=TRUE)
#' }
#'
#' @export paraline3D
paraline3D <- function(p,a,b,t)
{
  pname <-deparse(substitute(p))
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))

  if (!is.point(p,3) || !is.point(a,3) || !is.point(b,3))
  {stop('p, a and b must be numeric vectors of dimension 3')}

  if (!is.point(t,length(t)))
  {stop('t must be a numeric vector')}

  if (isTRUE(all.equal(a,b)))
  {stop('The lines are not well defined, the two points to define the line are concurrent!')}

  v<-b-a
  x0<-p[1]; y0<-p[2]; z0<-p[3];
  a<-v[1]; b<-v[2]; c<-v[3];

  x=x0+a*t
  y=y0+b*t
  z=z0+c*t

  vecs<-rbind(p,v)
  row.names(vecs) <- c("initial point","direction vector")

  line.desc<-paste("the line passing through point", pname, "parallel to the line joining points", aname, "and", bname)
  main.txt<-paste("Line Crossing point", pname, "\n in the direction of", bname, "-", aname)

  conda= c(a==0,a==1,a==-1,!(a==0|abs(a)==1))
  condb= c(b==0,b==1,b==-1,!(b==0|abs(b)==1))
  condc= c(c==0,c==1,c==-1,!(c==0|abs(c)==1))
  cx =switch(which(conda==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(a,"t",sep=""))
  cy =switch(which(condb==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(b,"t",sep=""))
  cz =switch(which(condc==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(c,"t",sep=""))
  cx0 =ifelse(x0==0,"",x0)
  cy0 =ifelse(y0==0,"",y0)
  cz0 =ifelse(z0==0,"",z0)
  signa=ifelse(sign(a)<=0,""," + ")
  signb=ifelse(sign(b)<=0,""," + ")
  signc=ifelse(sign(c)<=0,""," + ")

  eqn<- rbind(paste("x = ",cx0,signa,cx,sep=""),
              paste("y = ",cy0,signb,cy,sep=""),
              paste("z = ",cz0,signc,cz,sep=""))

  res<-list(
    desc=line.desc,
    mtitle=main.txt,
    pts=rbind(p,a,b),
    pnames=c(pname,aname,bname),
    vecs=vecs,
    vec.names=c(pname,paste(bname,"-",aname)),
    x=x,
    y=y,
    z=z,
    tsq=t,
    equation=eqn
  )
  class(res)<-"Lines3D"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The line crossing the 3D point \code{p} and perpendicular to the plane spanned by 3D points \code{a},
#' \code{b}, and \code{c}
#'
#' @description
#' An object of class \code{"Lines3D"}.
#' Returns the equation, \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinates of the line crossing 3D point \code{p} and perpendicular to the plane
#' spanned by 3D points \code{a}, \code{b}, and \code{c} (i.e., the line is in the direction of normal \code{vector} of this plane)
#' with the parameter \code{t} being provided in \code{vector} \code{t}.
#'
#' @param p A 3D point through which the straight line passes.
#' @param a,b,c 3D points which determine the plane to which the line passing through point \code{p} would be
#' perpendicular (i.e., the normal \code{vector} of this plane determines the direction of the straight line
#' passing through \code{p}).
#' @param t A scalar or a \code{vector} of scalars representing the parameter of the coordinates of the line
#' (for the form: \eqn{x=p_0 + A t}, \eqn{y=y_0 + B t}, and \eqn{z=z_0 + C t} where \eqn{p=(p_0,y_0,z_0)}
#' and normal vector\eqn{=(A,B,C)}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A description of the line}
#' \item{mtitle}{The \code{"main"} title for the plot of the line}
#' \item{points}{The input points that determine the line and plane, line crosses point \code{p} and plane is determined
#' by 3D points \code{a}, \code{b}, and \code{c}.}
#' \item{pnames}{The names of the input points that determine the line and plane; line would be perpendicular
#' to the plane.}
#' \item{vecs}{The point \code{p} and normal vector.}
#' \item{vec.names}{The names of the point \code{p} and the second entry is "normal vector".}
#' \item{x,y,z}{The \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinates of the point(s) of interest on the line perpendicular to the plane
#' determined by points \code{a}, \code{b}, and \code{c}.}
#' \item{tsq}{The scalar or the \code{vector} of the parameter in defining each coordinate of the line for the form:
#' \eqn{x=p_0 + A t}, \eqn{y=y_0 + B t}, and \eqn{z=z_0 + C t} where \eqn{p=(p_0,y_0,z_0)} and normal vector\eqn{=(A,B,C)}.}
#' \item{equation}{Equation of the line passing through point \code{p} and perpendicular to the plane determined by
#' points \code{a}, \code{b}, and \code{c} (i.e., line is in the direction of the normal \code{vector} N of the plane). The line equation
#' is in the form: \eqn{x=p_0 + A t}, \eqn{y=y_0 + B t}, and \eqn{z=z_0 + C t} where \eqn{p=(p_0,y_0,z_0)} and normal vector\eqn{=(A,B,C)}.}
#'
#' @seealso \code{\link{Line3D}}, \code{\link{paraline3D}} and \code{\link{perpline}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P<-c(1,1,1); Q<-c(1,10,4); R<-c(1,1,3); S<-c(3,9,12)
#'
#' cf<-as.numeric(Plane(Q,R,S,1,1)$coeff)
#' a<-cf[1]; b<-cf[2]; c<- -1;
#'
#' vecs<-rbind(Q,c(a,b,c))
#' pts<-rbind(P,Q,R,S)
#' perpline2plane(P,Q,R,S,.1)
#'
#' tr<-range(pts,vecs);
#' tf<-(tr[2]-tr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' tsq<-seq(-tf*10-tf,tf*10+tf,l=5)  #try also l=10, 20, or 100
#'
#' pln2pl<-perpline2plane(P,Q,R,S,tsq)
#' pln2pl
#' summary(pln2pl)
#' plot(pln2pl,theta = 225, phi = 30, expand = 0.7, facets = FALSE, scale = TRUE)
#'
#' xc<-pln2pl$x
#' yc<-pln2pl$y
#' zc<-pln2pl$z
#'
#' zr<-range(zc)
#' zf<-(zr[2]-zr[1])*.2
#' Rv<- -c(a,b,c)*zf*5
#'
#' Dr<-(Q+R+S)/3
#'
#' pts2<-rbind(Q,R,S)
#' xr<-range(pts2[,1],xc); yr<-range(pts2[,2],yc)
#' xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' yf<-(yr[2]-yr[1])*.1 #how far to go at the lower and upper ends in the y-coordinate
#' xs<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#' ys<-seq(yr[1]-yf,yr[2]+yf,l=5)  #try also l=10, 20, or 100
#'
#' plQRS<-Plane(Q,R,S,xs,ys)
#' z.grid<-plQRS$z
#'
#' Xlim<-range(xc,xs,pts[,1])
#' Ylim<-range(yc,ys,pts[,2])
#' Zlim<-range(zc,z.grid,pts[,3])
#'
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::persp3D(z = z.grid, x = xs, y = ys, theta =225, phi = 30,
#' main="Line Crossing P and \n Perpendicular to the Plane Defined by Q, R, S",
#' col="lightblue", ticktype = "detailed",
#'         xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.05,.05))
#'         #plane spanned by points Q, R, S
#' plot3D::lines3D(xc, yc, zc, bty = "g",pch = 20, cex = 2,col="red",
#' ticktype = "detailed",add=TRUE)
#' plot3D::arrows3D(Dr[1],Dr[2],Dr[3],Dr[1]+Rv[1],Dr[2]+Rv[2],Dr[3]+Rv[3], add=TRUE)
#' plot3D::points3D(pts[,1],pts[,2],pts[,3],add=TRUE)
#' plot3D::text3D(pts[,1],pts[,2],pts[,3],labels=c("P","Q","R","S"),add=TRUE)
#' plot3D::arrows3D(P[1],P[2],P[3]-zf,P[1],P[2],P[3],lty=2, add=TRUE)
#' plot3D::text3D(P[1],P[2],P[3]-zf,labels="initial point",add=TRUE)
#' plot3D::text3D(P[1],P[2],P[3]+zf/2,labels="P",add=TRUE)
#' plot3D::arrows3D(Dr[1],Dr[2],Dr[3],Dr[1]+Rv[1]/2,Dr[2]+Rv[2]/2,Dr[3]+Rv[3]/2,lty=2, add=TRUE)
#' plot3D::text3D(Dr[1]+Rv[1]/2,Dr[2]+Rv[2]/2,Dr[3]+Rv[3]/2,labels="normal vector",add=TRUE)
#' }
#' @export perpline2plane
perpline2plane <- function(p,a,b,c,t)
{
  pname <-deparse(substitute(p))
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))
  cname <-deparse(substitute(c))

  if (!is.point(p,3) || !is.point(a,3) || !is.point(b,3) || !is.point(c,3))
  {stop('p, a, b, and c must be numeric vectors of dimension 3')}

  if (!is.point(t,length(t)))
  {stop('t must be a numeric vector')}

  dab<-Dist(a,b); dac<-Dist(a,c); dbc<-Dist(b,c);
  sd<-sort(c(dab,dac,dbc))

  if (isTRUE(all.equal(sd[3] , sd[1]+sd[2])))
  {stop('The three points, a, b, and c are collinear; plane is not well-defined')}

  cf<-as.numeric(Plane(a,b,c,1,1)$coeff)
  A<-cf[1]; B<-cf[2]; C<- -1;
  x0<-p[1]; y0<-p[2]; z0<-p[3];

  x=x0+A*t
  y=y0+B*t
  z=z0+C*t

  vecs<-rbind(p,c(A,B,C))
  row.names(vecs) <- c("initial point","normal vector")

  line.desc<-paste("the line crossing point ", pname, " perpendicular to the plane spanned by points ", aname, ", ", bname, "and ", cname,sep="")
  main.txt<-paste("Line Crossing Point ", pname, "\n Perpendicular to the Plane\n Spanned by Points ", aname, ", ", bname, " and ", cname,sep="")

  condA= c(A==0,A==1,A==-1,!(A==0|abs(A)==1))
  condB= c(B==0,B==1,B==-1,!(B==0|abs(B)==1))
  condC= c(C==0,C==1,C==-1,!(C==0|abs(C)==1))
  cx =switch(which(condA==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(A,"t",sep=""))
  cy =switch(which(condB==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(B,"t",sep=""))
  cz =switch(which(condC==TRUE),"",paste("t",sep=""),
             paste("-t",sep=""),paste(C,"t",sep=""))
  cx0 =ifelse(x0==0,"",x0)
  cy0 =ifelse(y0==0,"",y0)
  cz0 =ifelse(z0==0,"",z0)
  signA=ifelse(sign(A)<=0,""," + ")
  signB=ifelse(sign(B)<=0,""," + ")
  signC=ifelse(sign(C)<=0,""," + ")

  eqn<- rbind(paste("x = ",cx0,signA,cx,sep=""),
              paste("y = ",cy0,signB,cy,sep=""),
              paste("z = ",cz0,signC,cz,sep=""))

  res<-list(
    desc=line.desc,
    mtitle=main.txt,
    pts=rbind(p,a,b,c),
    pnames=c(pname,aname,bname,cname),
    vecs=vecs,
    vec.names=c(pname,"normal vector"),
    x=x,
    y=y,
    z=z,
    tsq=t,
    equation=eqn
  )
  class(res)<-"Lines3D"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The point of intersection of a line and a plane
#'
#' @description Returns the point of the intersection of the line determined by the 3D points \eqn{p_1} and \eqn{p_2} and the plane spanned
#' by 3D points \code{p3}, \code{p4}, and \code{p5}.
#'
#' @param p1,p2 3D points that determine the straight line (i.e., through which the straight line passes).
#' @param p3,p4,p5 3D points that determine the plane (i.e., through which the plane passes).
#'
#' @return The coordinates of the point of intersection the line determined by the 3D points \eqn{p_1} and \eqn{p_2} and the
#' plane determined by 3D points \code{p3}, \code{p4}, and \code{p5}.
#'
#' @seealso \code{\link{int.2lines}} and \code{\link{int.circ.line}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' L1<-c(2,4,6); L2<-c(1,3,5);
#' A<-c(1,10,3); B<-c(1,1,3); C<-c(3,9,12)
#'
#' Pint<-int.line.plane(L1,L2,A,B,C)
#' Pint
#' pts<-rbind(L1,L2,A,B,C,Pint)
#'
#' tr<-max(Dist(L1,L2),Dist(L1,Pint),Dist(L2,Pint))
#' tf<-tr*1.1 #how far to go at the lower and upper ends in the x-coordinate
#' tsq<-seq(-tf,tf,l=5)  #try also l=10, 20, or 100
#'
#' lnAB3D<-Line3D(L1,L2,tsq)
#' xl<-lnAB3D$x
#' yl<-lnAB3D$y
#' zl<-lnAB3D$z
#'
#' xr<-range(pts[,1]); yr<-range(pts[,2])
#' xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' yf<-(yr[2]-yr[1])*.1 #how far to go at the lower and upper ends in the y-coordinate
#' xp<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#' yp<-seq(yr[1]-yf,yr[2]+yf,l=5)  #try also l=10, 20, or 100
#'
#' plABC<-Plane(A,B,C,xp,yp)
#' z.grid<-plABC$z
#'
#' res<-persp(xp,yp,z.grid, xlab="x",ylab="y",zlab="z",theta = -30, phi = 30, expand = 0.5,
#' col = "lightblue", ltheta = 120, shade = 0.05, ticktype = "detailed")
#' lines (trans3d(xl, yl, zl, pmat = res), col = 3)
#'
#' Xlim<-range(xl,pts[,1])
#' Ylim<-range(yl,pts[,2])
#' Zlim<-range(zl,pts[,3])
#'
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::persp3D(z = z.grid, x = xp, y = yp, theta =225, phi = 30, ticktype = "detailed"
#' ,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.1,.1),
#' expand = 0.7, facets = FALSE, scale = TRUE)
#'         #plane spanned by points A, B, C
#' #add the defining points
#' plot3D::points3D(pts[,1],pts[,2],pts[,3], pch = ".", col = "black", bty = "f", cex = 5,add=TRUE)
#' plot3D::points3D(Pint[1],Pint[2],Pint[3], pch = "*", col = "red", bty = "f", cex = 5,add=TRUE)
#' plot3D::lines3D(xl, yl, zl, bty = "g", cex = 2, ticktype = "detailed",add=TRUE)
#' }
#'
#' @export int.line.plane
int.line.plane <- function(p1,p2,p3,p4,p5)
{
  if (!is.point(p1,3) || !is.point(p2,3) || !is.point(p3,3) || !is.point(p4,3) || !is.point(p5,3))
  {stop('all arguments must be numeric vectors/points of dimension 3')}

  if (isTRUE(all.equal(p1,p2)))
  {stop('The line is not well defined.')}

  d34<-Dist(p3,p4); d35<-Dist(p3,p5); d45<-Dist(p4,p5);
  sd<-sort(c(d34,d35,d45))

  if (isTRUE(all.equal(sd[3] , sd[1]+sd[2])))
  {stop('the three points, p3, p4 and p5 are collinear; plane is not well-defined.')}

  num<- rbind(rep(1,4),cbind(p3,p4,p5,p1))
  denom<-rbind(c(1,1,1,0),cbind(p3,p4,p5,p2-p1))
  t<- - det(num)/det(denom)

  if (abs(t)==Inf)
  {stop('The line and plane are parallel, hence they do not intersect.')}

  v<-p2-p1
  x0<-p1[1]; y0<-p1[2]; z0<-p1[3];
  a<-v[1]; b<-v[2]; c<-v[3];

  x=x0+a*t
  y=y0+b*t
  z=z0+c*t

  c(x,y,z)
} #end of the function
#'

#################################################################

#' @title The plane passing through three distinct 3D points \code{a}, \code{b}, and \code{c}
#'
#' @description
#' An object of class \code{"Planes"}.
#' Returns the equation and \eqn{z}-coordinates of the plane passing through three distinct 3D points \code{a}, \code{b}, and \code{c}
#' with \eqn{x}- and \eqn{y}-coordinates are provided in vectors \code{x} and \code{y}, respectively.
#'
#' @param a,b,c 3D points that determine the plane (i.e., through which the plane is passing).
#' @param x,y Scalars or vectors of scalars representing the \eqn{x}- and \eqn{y}-coordinates of the plane.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A description of the plane}
#' \item{points}{The input points \code{a}, \code{b}, and \code{c} through which the plane is passing
#' (stacked row-wise, i.e., row 1 is point \code{a}, row 2 is point \code{b} and row 3 is point \code{c}).}
#' \item{x,y}{The input vectors which constitutes the \eqn{x}- and \eqn{y}-coordinates of the point(s) of interest on the
#' plane. \code{x} and \code{y} can be scalars or vectors of scalars.}
#' \item{z}{The output \code{vector} which constitutes the \eqn{z}-coordinates of the point(s) of interest on the plane.
#' If \code{x} and \code{y} are scalars, \code{z} will be a scalar and
#' if \code{x} and \code{y} are vectors of scalars, then \code{z} needs to be a \code{matrix} of scalars,
#' containing the \eqn{z}-coordinate for each pair of \code{x} and \code{y} values.}
#' \item{coeff}{Coefficients of the plane (in the \eqn{z = A x+B y+C} form).}
#' \item{equation}{Equation of the plane in long form}
#' \item{equation2}{Equation of the plane in short form, to be inserted on the plot}
#'
#' @seealso \code{\link{paraplane}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P1<-c(1,10,3); P2<-c(1,1,3); P3<-c(3,9,12) #also try P2=c(2,2,3)
#'
#' pts<-rbind(P1,P2,P3)
#' Plane(P1,P2,P3,.1,.2)
#'
#' xr<-range(pts[,1]); yr<-range(pts[,2])
#' xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' yf<-(yr[2]-yr[1])*.1 #how far to go at the lower and upper ends in the y-coordinate
#' x<-seq(xr[1]-xf,xr[2]+xf,l=5)  #try also l=10, 20, or 100
#' y<-seq(yr[1]-yf,yr[2]+yf,l=5)  #try also l=10, 20, or 100
#'
#' plP123<-Plane(P1,P2,P3,x,y)
#' plP123
#' summary(plP123)
#' plot(plP123,theta = 225, phi = 30, expand = 0.7, facets = FALSE, scale = TRUE)
#'
#' z.grid<-plP123$z
#'
#' persp(x,y,z.grid, xlab="x",ylab="y",zlab="z",
#' theta = -30, phi = 30, expand = 0.5, col = "lightblue",
#'       ltheta = 120, shade = 0.05, ticktype = "detailed")
#'
#' zr<-max(z.grid)-min(z.grid)
#' Pts<-rbind(P1,P2,P3)+rbind(c(0,0,zr*.1),c(0,0,zr*.1),c(0,0,zr*.1))
#' Mn.pts<-apply(Pts,2,mean)
#'
#' plot3D::persp3D(z = z.grid, x = x, y = y,theta = 225, phi = 30, expand = 0.3,
#' main = "Plane Crossing Points P1, P2, and P3", facets = FALSE, scale = TRUE)
#' #plane spanned by points P1, P2, P3
#' #add the defining points
#' plot3D::points3D(Pts[,1],Pts[,2],Pts[,3], add=TRUE)
#' plot3D::text3D(Pts[,1],Pts[,2],Pts[,3], c("P1","P2","P3"),add=TRUE)
#' plot3D::text3D(Mn.pts[1],Mn.pts[2],Mn.pts[3],plP123$equation,add=TRUE)
#' #plot3D::polygon3D(Pts[,1],Pts[,2],Pts[,3], add=TRUE)
#' }
#'
#' @export Plane
Plane <- function(a,b,c,x,y)
{
  aname <-deparse(substitute(a))
  bname <-deparse(substitute(b))
  cname <-deparse(substitute(c))

  if (!is.point(a,3) || !is.point(b,3) || !is.point(c,3))
  {stop('a, b,and c must be numeric points of dimension 3')}

  if (!is.point(x,length(x)) || !is.point(y,length(y)))
  {stop('x and y must be numeric vectors')}

  dab<-Dist(a,b); dac<-Dist(a,c); dbc<-Dist(b,c);
  sd<-sort(c(dab,dac,dbc))

  if (isTRUE(all.equal(sd[3] , sd[1]+sd[2])))
  {stop('the three points, a, b, and c, are collinear; plane is not well-defined')}

  mat<-rbind(a,b,c)
  D<-det(mat)
  if (D==0)
  {stop('The three points are collinear (i.e., they lie on a plane or concur at the same point')}

  matx<-cbind(rep(1,3),mat[,2],mat[,3])
  maty<-cbind(mat[,1],rep(1,3),mat[,3])
  matz<-cbind(mat[,1],mat[,2],rep(1,3))
  a0<-(-1/D)*det(matx)
  b0<-(-1/D)*det(maty)
  c0<-(-1/D)*det(matz)

  A<- -a0/c0; B<- -b0/c0; C<- -1/c0;

  z.grid <- outer(x, y, function(a, b) A*a+B*b+C)

  plane.desc<-paste("Plane Passing through Points ", aname,", ", bname,",", " and ", cname,sep="")
  plane.title<-paste("Plane Crossing Points ", aname,", ", bname,",", " and ", cname,sep="")

  pts<-rbind(a,b,c)
  row.names(pts)<-c(aname,bname,cname)

  coeffs<-c(A,B,C)
  names(coeffs)<-c("A","B","C")

  condA= c(A==0,A==1,!(A==0|A==1))
  condB= c(B==0,B==1,!(B==0|B==1))
  compx =switch(which(condA==TRUE),"",paste("x",sep=""),paste(A,"x",sep=""))
  compy =switch(which(condB==TRUE),"",paste("y",sep=""),paste(B,"y",sep=""))

  signB=ifelse(sign(B)<=0,""," + ")
  signC=ifelse(sign(C)<=0," "," + ")

  eqn=paste("z = ",compx,signB,compy,signC,C,sep="");

  compp2 =switch(which(condA==TRUE),"",paste("x",sep=""),paste(my_round(A,2),"x",sep=""))
  compy2 =switch(which(condB==TRUE),"",paste("y",sep=""),paste(my_round(B,2),"y",sep=""))
  eqn2=paste("z = ",compp2,signB,compy2,signC,my_round(C,2),sep="");

  res<-list(
    desc=plane.desc,
    main.title=plane.title,
    points=pts,
    x=x,
    y=y,
    z=z.grid,
    coeff=coeffs,
    equation=eqn,
    equation2=eqn2
    )

  class(res)<-"Planes"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The distance from a point to a plane spanned by three 3D points
#'
#' @description Returns the distance from a point \code{p} to the plane passing through points \code{a}, \code{b}, and \code{c} in 3D space.
#'
#' @param p A 3D point, distance from \code{p} to the plane passing through points \code{a}, \code{b},
#' and \code{c} are to be computed.
#' @param a,b,c 3D points that determine the plane (i.e., through which the plane is passing).
#'
#' @return A \code{list} with two elements
#' \item{dis}{Distance from point \code{p} to the plane spanned by 3D points \code{a}, \code{b}, and \code{c}}
#' \item{cl2pl}{The closest point on the plane spanned by 3D points \code{a}, \code{b}, and \code{c} to the point \code{p}}
#'
#' @seealso \code{\link{dist.point2line}}, \code{\link{dist.point2set}}  and \code{\link{Dist}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P<-c(5,2,40)
#' P1<-c(1,2,3); P2<-c(3,9,12); P3<-c(1,1,3);
#'
#' dis<-dist.point2plane(P,P1,P2,P3);
#' dis
#' Pr<-dis$proj #projection on the plane
#'
#' xseq<-seq(0,10,l=5)  #try also l=10, 20, or 100
#' yseq<-seq(0,10,l=5)  #try also l=10, 20, or 100
#'
#' pl.grid<-Plane(P1,P2,P3,xseq,yseq)$z
#'
#' plot3D::persp3D(z = pl.grid, x = xseq, y = yseq, theta =225, phi = 30, ticktype = "detailed",
#' expand = 0.7, facets = FALSE, scale = TRUE,
#' main="Point P and its Orthogonal Projection \n on the Plane Defined by P1, P2, P3")
#' #plane spanned by points P1, P2, P3
#' #add the vertices of the tetrahedron
#' plot3D::points3D(P[1],P[2],P[3], add=TRUE)
#' plot3D::points3D(Pr[1],Pr[2],Pr[3], add=TRUE)
#' plot3D::segments3D(P[1], P[2], P[3], Pr[1], Pr[2],Pr[3], add=TRUE,lwd=2)
#'
#' plot3D::text3D(P[1]-.5,P[2],P[3]+1, c("P"),add=TRUE)
#' plot3D::text3D(Pr[1]-.5,Pr[2],Pr[3]+2, c("Pr"),add=TRUE)
#'
#' persp(xseq,yseq,pl.grid, xlab="x",ylab="y",zlab="z",theta = -30,
#' phi = 30, expand = 0.5, col = "lightblue",
#'       ltheta = 120, shade = 0.05, ticktype = "detailed")
#'
#' }
#'
#' @export dist.point2plane
dist.point2plane <- function(p,a,b,c)
{
  if (!is.point(p,3) || !is.point(a,3) || !is.point(b,3) || !is.point(c,3))
  {stop('all arguments must be numeric points of dimension 3')}

  p1<-p[1]; p2<-p[2]; p3<-p[3];
  a1<-a[1]; a2<-a[2]; a3<-a[3];
  b1<-b[1]; b2<-b[2]; b3<-b[3];
  c1<-c[1]; c2<-c[2]; c3<-c[3];

  Cx <-(-a2*b3+a2*c3+a3*b2-a3*c2-b2*c3+b3*c2)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1);
  Cy <-(a1*b3-a1*c3-a3*b1+a3*c1+b1*c3-b3*c1)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1);
  Cz <--1;
  if (abs(Cx)==Inf || abs(Cy)==Inf)  #this part is for when the plane is vertical so Cx or Cy are Inf or -Inf
  {ra<-range(a); a<-a+(ra[2]-ra[1])*runif(3,0,10^(-10))
  rb<-range(b); b<-b+(rb[2]-rb[1])*runif(3,0,10^(-10))
  rc<-range(c); c<-c+(rc[2]-rc[1])*runif(3,0,10^(-10))
  a1<-a[1]; a2<-a[2]; a3<-a[3];
  b1<-b[1]; b2<-b[2]; b3<-b[3];
  c1<-c[1]; c2<-c[2]; c3<-c[3];
  Cx <-(-a2*b3+a2*c3+a3*b2-a3*c2-b2*c3+b3*c2)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1);
  Cy <-(a1*b3-a1*c3-a3*b1+a3*c1+b1*c3-b3*c1)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1);
  }
  C0 <-(a1*b2*c3-a1*b3*c2-a2*b1*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1)/(a1*b2-a1*c2-a2*b1+a2*c1+b1*c2-b2*c1);

  num.dis<-abs(Cx*p1+Cy*p2+Cz*p3+C0); #numerator of the distance equation

  t<-(Cx*a1-Cx*p1+Cy*a2-Cy*p2+Cz*a3-Cz*p3)/(Cx^2+Cy^2+Cz^2)
  prj<-p+t*c(Cx,Cy,Cz)

  dis<-num.dis/sqrt(Cx^2+Cy^2+Cz^2);

  list(distance=dis, #distance
       proj.pt2plane=prj #c(pr.x,pr.y,pr.z) #point of orthogonal projection (i.e. closest point) on the plane
  )
} #end of the function
#'
