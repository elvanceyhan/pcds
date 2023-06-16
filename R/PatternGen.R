#PatternGen.R;
#Contains functions to generate spatial point patterns in 2D regions
#################################################################

#' @title Generation of points segregated (in a radial or circular fashion)
#' from a given set of points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} 2D points uniformly
#' in \eqn{(a_1-e,a_1+e) \times (a_1-e,a_1+e) \setminus B(y_i,e)}
#' (\eqn{a_1} and \eqn{b1} are denoted as \code{a1} and \code{b1}
#' as arguments) where
#' \eqn{Y_p=(y_1,y_2,\ldots,y_{n_y})}
#' with \eqn{n_y} being number of \code{Yp} points
#' for various values of \code{e}
#' under the segregation pattern
#' and \eqn{B(y_i,e)} is the ball centered at \eqn{y_i}
#' with radius \code{e}.
#'
#' Positive values of \code{e} yield realizations from the segregation pattern
#' and nonpositive values of \code{e} provide
#' a type of complete spatial randomness (CSR),
#' \code{e} should not be too large
#' to make the support of generated points empty,
#' \code{a1} is defaulted
#' to the minimum of the \eqn{x}-coordinates of the \code{Yp} points,
#' \code{a2} is defaulted
#' to the maximum of the \eqn{x}-coordinates of the \code{Yp} points,
#' \code{b1} is defaulted
#' to the minimum of the \eqn{y}-coordinates of the \code{Yp} points,
#' \code{b2} is defaulted
#' to the maximum of the \eqn{y}-coordinates of the \code{Yp} points.
#'
#' @param n A positive integer representing the number of points to be generated.
#' @param Yp A set of 2D points representing the reference points.
#' The generated points are segregated
#' (in a circular or radial fashion) from these points.
#' @param e A positive real number
#' representing the radius of the balls centered at \code{Yp} points.
#' These balls are forbidden for the generated points
#' (i.e., generated points would be in the complement of union of these
#' balls).
#' @param a1,a2 Real numbers representing
#' the range of \eqn{x}-coordinates in the region
#' (default is the range of \eqn{x}-coordinates of the \code{Yp} points).
#' @param b1,b2 Real numbers representing
#' the range of \eqn{y}-coordinates in the region
#' (default is the range of \eqn{y}-coordinates of the \code{Yp} points).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{Radial
#' (i.e., circular) exclusion parameter of the segregation pattern}
#' \item{ref.points}{The input set of reference points \code{Yp},
#' i.e., points from which generated points are
#' segregated.}
#' \item{gen.points}{The output set of generated points segregated
#' from \code{Yp} points}
#' \item{tri.Yp}{Logical output for triangulation
#' based on \code{Yp} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Yp} points is
#' to be implemented (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and
#' the number of reference (i.e., \code{Yp}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the generated points}
#'
#' @seealso \code{\link{rassoc.circular}}, \code{\link{rseg.std.tri}},
#' \code{\link{rsegII.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-100; ny<-4;  #try also nx<-1000; ny<-10
#' e<-.15;  #try also e<- -.1; #a negative e provides a CSR realization
#' #with default bounding box (i.e., unit square)
#' Y<-cbind(runif(ny),runif(ny))
#'
#' Xdt<-rseg.circular(nx,Y,e)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' #with default bounding box (i.e., unit square)
#' Y<-cbind(runif(ny),runif(ny))
#' Xdt<-Xdt$gen.points
#' Xlim<-range(Xdt[,1],Y[,1]);
#' Ylim<-range(Xdt[,2],Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Y,asp=1,pch=16,col=2,lwd=2, xlab="x",ylab="y",
#'      main="Circular Segregation of X points from Y Points",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' points(Xdt)
#'
#' #with a rectangular bounding box
#' a1<-0; a2<-10;
#' b1<-0; b2<-5;
#' e<-1.5;
#' Y<-cbind(runif(ny,a1,a2),runif(ny,b1,b2))
#' #try also Y<-cbind(runif(ny,a1,a2/2),runif(ny,b1,b2/2))
#'
#' Xdt<-rseg.circular(nx,Y,e,a1,a2,b1,b2)$gen.points
#' Xlim<-range(Xdt[,1],Y[,1]); Ylim<-range(Xdt[,2],Y[,2])
#'
#' plot(Y,pch=16,asp=1,col=2,lwd=2, xlab="x",ylab="y",
#'      main="Circular Segregation of X points from Y Points",
#'      xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(Xdt)
#' }
#'
#' @export rseg.circular
rseg.circular <- function(n,Yp,e,a1=min(Yp[,1]),a2=max(Yp[,1]),b1=min(Yp[,2]),b2=max(Yp[,2]))
{
  if (is.point(Yp))
  { Yp<-matrix(Yp,ncol=2)
  } else
  {Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 )
  {stop('Yp must be of dimension nx2')}
  }

  ny<-nrow(Yp)

  if (!is.point(e,1))
  {stop('e must be a scalar')}

  if (ny>1 && e>=min((a2-a1)/2, (b2-b1)/2))
  {warning("warning: e is taken large, so it might take a very long time!")}

  if (e<=0)
  {warning("e is taken nonpositive, so the pattern is actually uniform in the study region!")}

  i<-0
  Xdt<-matrix(0,n,2)
  while (i<=n)
  {
    x<-c(runif(1,a1-e,a2+e),runif(1,b1-e,b2+e))
    if (dist.point2set(x,Yp)$distance>=e)
    {Xdt[i,]<-x;
    i<-i+1}
  }

  param<-e
  names(param)<-"exclusion parameter"
  typ<-paste("Segregation of ",n, " X points from ", ny," Y points with circular exclusion parameter e = ",e,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Segregation of Class X from Class Y"
  main.txt<-paste("Segregation of Two Classes \n with circular exclusion parameter e = ",e,sep="")
  res<-list(
    type=typ,
    parameters=param,
    gen.points=Xdt, #generated points segregated from Y points
    ref.points=Yp,
    #reference points, i.e., points from which generated points are segregated
    tri.Yp=FALSE,
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=c(a1-e,a2+e),
    ylimit=c(b1-e,b2+e)
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points associated (in a radial or circular fashion)
#' with a given set of points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} 2D points uniformly
#' in \eqn{(a_1-e,a_1+e) \times (a_1-e,a_1+e) \cap U_i B(y_i,e)}
#' (\eqn{a_1} and \eqn{b1} are denoted as
#' \code{a1} and \code{b1} as arguments)
#' where \eqn{Y_p=(y_1,y_2,\ldots,y_{n_y})} with
#' \eqn{n_y} being number of \code{Yp} points
#' for various values of \code{e} under the association pattern
#' and \eqn{B(y_i,e)} is the ball centered at \eqn{y_i} with radius \code{e}.
#'
#' \code{e} must be positive
#' and very large values of \code{e} provide patterns close to CSR.
#' \code{a1} is defaulted
#' to the minimum of the \eqn{x}-coordinates of the \code{Yp} points,
#' \code{a2} is defaulted
#' to the maximum of the \eqn{x}-coordinates of the \code{Yp} points,
#' \code{b1} is defaulted
#' to the minimum of the \eqn{y}-coordinates of the \code{Yp} points,
#' \code{b2} is defaulted
#' to the maximum of the \eqn{y}-coordinates of the \code{Yp} points.
#' This function is also very similar to \code{\link{rassoc.matern}},
#' where \code{rassoc.circular}
#' needs the study window to be specified,
#' while \code{\link{rassoc.matern}} does not.
#'
#' @param n A positive integer representing the number of points
#' to be generated.
#' @param Yp A set of 2D points representing the reference points.
#' The generated points are associated
#' (in a circular or radial fashion) with these points.
#' @param e A positive real number
#' representing the radius of the balls centered at \code{Yp} points.
#' Only these balls are allowed for the generated points
#' (i.e., generated points would be in the union of these balls).
#' @param a1,a2 Real numbers
#' representing the range of \eqn{x}-coordinates in the region
#' (default is the range of \eqn{x}-coordinates of the \code{Yp} points).
#' @param b1,b2 Real numbers
#' representing the range of \eqn{y}-coordinates in the region
#' (default is the range of \eqn{y}-coordinates of the \code{Yp} points).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{Radial attraction parameter of the association pattern}
#' \item{ref.points}{The input set of attraction points \code{Yp},
#' i.e., points with which generated points are associated.}
#' \item{gen.points}{The output set of generated points
#' associated with \code{Yp} points}
#' \item{tri.Yp}{Logical output for triangulation
#' based on \code{Yp} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Yp} points is
#' to be implemented (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number of attraction
#' (i.e., \code{Yp}) points.}
#' \item{xlimit,ylimit}{The possible range of the \eqn{x}-
#' and \eqn{y}-coordinates of the generated points.}
#'
#' @seealso \code{\link{rseg.circular}}, \code{\link{rassoc.std.tri}},
#' \code{\link{rassocII.std.tri}}, \code{\link{rassoc.matern}},
#' and \code{\link{rassoc.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-100; ny<-4;  #try also nx<-1000; ny<-10;
#'
#' e<-.15;
#' #with default bounding box (i.e., unit square)
#' Y<-cbind(runif(ny),runif(ny))
#'
#' Xdt<-rassoc.circular(nx,Y,e)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xdt<-Xdt$gen.points
#' Xlim<-range(Xdt[,1],Y[,1]);
#' Ylim<-range(Xdt[,2],Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Y,asp=1,xlab="x",ylab="y",
#' main="Circular Association of X points with Y Points",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01),
#'      pch=16,col=2,lwd=2)
#' points(Xdt)
#'
#' #with default bounding box (i.e., unit square)
#' Xlim<-range(Xdt[,1],Y[,1]);
#' Ylim<-range(Xdt[,2],Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Y,asp=1,xlab="x",ylab="y",
#' main="Circular Association of X points with Y Points",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01),pch=16,
#'      col=2,lwd=2)
#' points(Xdt)
#'
#' #with a rectangular bounding box
#' a1<-0; a2<-10;
#' b1<-0; b2<-5;
#' e<-1.1;  #try also e<-5; #pattern very close to CSR!
#'
#' Y<-cbind(runif(ny,a1,a2),runif(ny,b1,b2))
#' #try also Y<-cbind(runif(ny,a1,a2/2),runif(ny,b1,b2/2))
#'
#' Xdt<-rassoc.circular(nx,Y,e,a1,a2,b1,b2)$gen.points
#' Xlim<-range(Xdt[,1],Y[,1]);
#' Ylim<-range(Xdt[,2],Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Y,asp=1,xlab="x",ylab="y",
#' main="Circular Association of X points with Y Points",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01),
#'      pch=16,col=2,lwd=2)
#' points(Xdt)
#' }
#'
#' @export rassoc.circular
rassoc.circular <- function(n,Yp,e,a1=min(Yp[,1]),a2=max(Yp[,1]),b1=min(Yp[,2]),b2=max(Yp[,2]))
{
  if (is.point(Yp))
  { Yp<-matrix(Yp,ncol=2)
  } else
  {Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 )
  {stop('Yp must be of dimension nx2')}
  }

  ny<-nrow(Yp)

  if (!is.point(e,1) || e<=0)
  {stop('e must be a positive scalar')}

  if (ny>1 && e>=max((a2-a1)/2, (b2-b1)/2))
  {warning("e is taken too large, the pattern will be very close to CSR!")}

  i<-0
  Xdt<-matrix(0,n,2)
  while (i<=n)
  {
    x<-c(runif(1,a1-e,a2+e),runif(1,b1-e,b2+e))
    if (dist.point2set(x,Yp)$distance<=e)
    {Xdt[i,]<-x;
    i<-i+1}
  }
  param<-e
  names(param)<-"attraction parameter"
  typ<-paste("Association of ",n, " points with ", ny," Y points with circular attraction parameter e = ",e,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Association of one Class with Class Y"
  main.txt<-paste("Association of Two Classes \n with circular attraction parameter e = ",e,sep="")
  res<-list(
    type=typ,
    parameters=param,
    gen.points=Xdt, #generated points associated with Y points
    ref.points=Yp, #attraction points, i.e., points to which generated points are associated
    tri.Yp=FALSE,
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=c(a1-e,a2+e),
    ylimit=c(b1-e,b2+e)
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points associated (in a Matern-like fashion)
#' to a given set of points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} 2D points uniformly in \eqn{\cup B(y_i,e)}
#' where \eqn{Y_p=(y_1,y_2,\ldots,y_{n_y})}
#' with \eqn{n_y} being number of \code{Yp} points
#' for various values of \code{e} under the association pattern
#' and \eqn{B(y_i,e)} is the ball centered
#' at \eqn{y_i} with radius \code{e}.
#'
#' The pattern resembles the Matern cluster pattern
#' (see \code{\link[spatstat.random]{rMatClust}} in the
#' \code{spatstat.random} package for further information
#' (\insertCite{baddeley:2005;textual}{pcds}).
#' \code{rMatClust(kappa, scale, mu, win)} in the simplest
#' case generates a uniform Poisson point process of "parent" points
#' with intensity \code{kappa}.
#' Then each parent point is replaced by a random cluster of
#' "offspring" points, the number of points per cluster
#' being Poisson(\code{mu}) distributed,
#' and their positions being placed
#' and uniformly inside a disc of radius scale centered on the parent point.
#' The resulting point pattern is a realization of the classical
#' "stationary Matern cluster process" generated inside the
#' window \code{win}.
#'
#' The main difference of \code{rassoc.matern}
#' and \code{\link[spatstat.random]{rMatClust}}
#' is that the parent points are \code{Yp} points
#' which are given beforehand
#' and we do not discard them in the end in \code{rassoc.matern}
#' and the offspring points are the points associated
#' with the reference points, \code{Yp};
#' \code{e} must be positive and very large values of \code{e}
#' provide patterns close to CSR.
#'
#' This function is also very similar to \code{\link{rassoc.circular}},
#' where \code{\link{rassoc.circular}} needs the study window to be specified,
#' while \code{rassoc.matern} does not.
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param Yp A set of 2D points representing the reference points.
#' The generated points are associated
#' (in a Matern-cluster like fashion) with these points.
#' @param e A positive real number representing the radius of the balls
#' centered at \code{Yp} points.
#' Only these balls are allowed for the generated points
#' (i.e., generated points would be in the union of these balls).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{Radial (i.e., circular) attraction parameter
#' of the association pattern.}
#' \item{ref.points}{The input set of attraction points \code{Yp},
#' i.e., points with which generated points are associated.}
#' \item{gen.points}{The output set of generated points associated
#' with \code{Yp} points.}
#' \item{tri.Yp}{Logical output for triangulation
#' based on \code{Yp} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Yp} points is
#' to be implemented (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number of
#' attraction (i.e., \code{Yp}) points.}
#' \item{xlimit,ylimit}{The possible ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the generated points.}
#'
#' @seealso \code{\link{rassoc.circular}}, \code{\link{rassoc.std.tri}},
#' \code{\link{rassocII.std.tri}}, \code{\link{rassoc.multi.tri}},
#' \code{\link{rseg.circular}}, and \code{\link[spatstat.random]{rMatClust}}
#' in the \code{spatstat.random} package
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-100; ny<-4;  #try also nx<-1000; ny<-10;
#'
#' e<-.15;
#'  #try also e<-1.1; #closer to CSR than association, as e is large
#'
#' #Y points uniform in unit square
#' Y<-cbind(runif(ny),runif(ny))
#'
#' Xdt<-rassoc.matern(nx,Y,e)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xdt<-Xdt$gen.points
#' Xlim<-range(Xdt[,1],Y[,1]);
#' Ylim<-range(Xdt[,2],Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Y,asp=1,xlab="x",ylab="y",
#' main="Matern-like Association of X points with Y Points",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01),
#'      pch=16,col=2,lwd=2)
#' points(Xdt)
#'
#' a1<-0; a2<-10;
#' b1<-0; b2<-5;
#' e<-1.1;
#'
#' #Y points uniform in a rectangle
#' Y<-cbind(runif(ny,a1,a2),runif(ny,b1,b2))
#' #try also Y<-cbind(runif(ny,a1,a2/2),runif(ny,b1,b2/2))
#'
#' Xdt<-rassoc.matern(nx,Y,e)$gen.points
#' Xlim<-range(Xdt[,1],Y[,1]);
#' Ylim<-range(Xdt[,2],Y[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Y,asp=1,xlab="x",ylab="y",
#' main="Matern-like Association of X points with Y Points",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01),pch=16,col=2,lwd=2)
#' points(Xdt)
#' }
#'
#' @export rassoc.matern
rassoc.matern <- function(n,Yp,e)
{
  if (is.point(Yp))
  { Yp<-matrix(Yp,ncol=2)
  } else
  {Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 )
  {stop('Yp must be of dimension nx2')}
  }

  if (!is.point(e,1) || e<=0)
  {stop('e must be a positive scalar')}

  ny<-nrow(Yp)
  indy<-sample(1:ny, n, replace = TRUE)
  Xdt<-matrix(0,n,2)

  for (i in 1:n)
  {
    Yp2<-Yp[indy[i],]
    t<-2*pi*runif(1)
    r<-e*sqrt(runif(1))
    Xdt[i,]<-Yp2 + r*c(cos(t),sin(t))
  }

  param<-e
  names(param)<-"attraction parameter"
  typ<-paste("Matern-like Association of ",n, " points with ", ny," Y points with circular attraction parameter e = ",e,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Matern-like Association of one Class with Class Y"
  main.txt<-paste("Matern-like Association of Two Classes \n Circular Attraction Parameter e = ",e,sep="")
  Xlim<-range(Yp[,1])+c(-e,e)
  Ylim<-range(Yp[,2])+c(-e,e)
  res<-list(
    type=typ,
    parameters=param,
    gen.points=Xdt, #generated points associated with Y points
    ref.points=Yp, #attraction points, i.e., points to which generated points are associated
    tri.Yp=FALSE,
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    num.points=npts,
    xlimit=Xlim,
    ylimit=Ylim
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Uniform Points in the first one-sixth of
#' standard equilateral triangle
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly
#' in the first 1/6th of the standard equilateral triangle \eqn{T_e=(A,B,C)}
#' with vertices with \eqn{A=(0,0)}; \eqn{B=(1,0)}, \eqn{C=(1/2,\sqrt{3}/2)}
#' (see the examples below).
#' The first 1/6th of the standard equilateral triangle is the triangle with vertices
#' \eqn{A=(0,0)}, \eqn{(1/2,0)}, \eqn{C=(1/2,\sqrt{3}/6)}.
#'
#' @param n a positive integer representing number of uniform points
#' to be generated
#' in the first one-sixth of \eqn{T_e}.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{support}{The vertices of the support of
#' the uniformly generated points}
#' \item{gen.points}{The output set of uniformly generated points
#' in the first 1/6th of
#' the standard equilateral triangle.}
#' \item{out.region}{The outer region for the one-sixth of \eqn{T_e},
#' which is just \eqn{T_e} here.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points
#' and the number of vertices of the support (i.e., \code{Y}) points.}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the generated,
#' support and outer region points}
#'
#' @seealso \code{\link{runif.std.tri}}, \code{\link{runif.basic.tri}},
#' \code{\link{runif.tri}}, and \code{\link{runif.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' CM<-(A+B+C)/3;
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' nx<-100  #try also nx<-1000
#'
#' #data generation step
#' set.seed(1)
#' Xdt<-runif.std.tri.onesixth(nx)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xd<-Xdt$gen.points
#'
#' #plot of the data with the regions in the equilateral triangle
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,asp=1,pch=".",xlim=Xlim+xd*c(-.01,.01),
#' ylim=Ylim+yd*c(-.01,.01),xlab=" ",ylab=" ",
#'      main="first 1/6th of the \n standard equilateral triangle")
#' polygon(Te)
#' L<-Te; R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' polygon(rbind(A,D3,CM),col=5)
#' points(Xd)
#'
#' #letter annotation of the plot
#' txt<-rbind(A,B,C,CM,D1,D2,D3)
#' xc<-txt[,1]+c(-.02,.02,.02,.04,.05,-.03,0)
#' yc<-txt[,2]+c(.02,.02,.02,.03,0,.03,-.03)
#' txt.str<-c("A","B","C","CM","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
runif.std.tri.onesixth <- function(n)
{
  X <-matrix(0,n,2);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,0,.5); y<-runif(1,0,0.2886751347);
    if (y<.5773502694*x)
    {X[i,]<-c(x,y);
    ct<-1;
    }
  }
  }
  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2)
  CM<-(A+B+C)/3
  Te<-rbind(A,B,C)  #std eq triangle
  Te1_6<-rbind(A,(A+B)/2,CM)  #std eq triangle

  typ<-"Uniform Distribution in one-sixth of the standard equilateral triangle"
  txt<-paste(n, " uniform points in one-sixth of the standard equilateral triangle")
  main.txt<-"Uniform points in one-sixth of the\n standard equilateral triangle"

  npts<-c(n,3)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], "is the number of Uniform points \n",names(npts)[2], "is the number of vertices of the support region")
  row.names(Te1_6)<-c()

  res<-list(
    type=typ,
    gen.points=X, #uniformly generated points in one-sixth of std eq triangle
    tess.points=Te1_6, #tessellation points whose convex hull constitutes the support of the uniform points
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    out.region=Te, #outer region for the one-sixth of Te
    num.points=npts,  txt4pnts=txt4pnts,
    xlimit=range(Te[,1]),
    ylimit=range(Te[,2])
  )

  class(res)<-"Uniform"

  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Generation of Uniform Points in the Standard Equilateral Triangle
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly
#' in the standard equilateral triangle \eqn{T_e=T(A,B,C)}
#' with vertices \eqn{A=(0,0)}, \eqn{B=(1,0)}, and \eqn{C=(1/2,\sqrt{3}/2)}.
#'
#' @param n A positive integer representing the number of uniform points
#' to be generated in the standard equilateral triangle \eqn{T_e}.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{tess.points}{The vertices of the support region of
#' the uniformly generated points, it is the
#' standard equilateral triangle \eqn{T_e} for this function}
#' \item{gen.points}{The output set of generated points uniformly
#' in the standard equilateral triangle \eqn{T_e}.}
#' \item{out.region}{The outer region which contains the support region,
#' \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern
#' from which points are to be generated}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number
#' of vertices of the support points (here it is 3).}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the support, \eqn{T_e}}
#'
#' @seealso \code{\link{runif.basic.tri}}, \code{\link{runif.tri}},
#' and \code{\link{runif.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' n<-100
#'
#' set.seed(1)
#' Xdt<-runif.std.tri(n)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' Xp<-runif.std.tri(n)$gen.points
#' plot(Te,asp=1,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.01,.01),
#' ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp)
#' }
#'
#' @export runif.std.tri
runif.std.tri <- function(n)
{
  Xdt <-matrix(0,n,2);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,0,1); y<-runif(1,0,.8660254040);
    if (y<1.732050808*x && y<1.732050808-1.732050808*x)
    {Xdt[i,]<-c(x,y);
    ct<-1;
    }
  }
  }
  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2)
  Te<-rbind(A,B,C)  #std eq triangle

  typ<-"Uniform Distribution in the Standard Equilateral Triangle"
  main.txt<-"Uniform Points in the Standard Equilateral Triangle"
  txt<-paste(n, " uniform points in the standard equilateral triangle")

  npts<-c(n,3)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], "is the number of uniform points \n",names(npts)[2], "is the number of vertices of the support region")

  res<-list(
    type=typ,
    mtitle=main.txt,
    tess.points=Te,
    #tessellation points whose convex hull constitutes the support of the uniform points
    gen.points=Xdt, #uniformly generated points in the std equilateral triangle
    out.region=NULL, #outer region for Te
    desc.pat=txt, #description of the pattern
    num.points=npts,  txt4pnts=txt4pnts,
    xlimit=range(Te[,1]),
    ylimit=range(Te[,2])
  )

  class(res)<-"Uniform"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Uniform Points in the standard basic triangle
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly
#' in the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0}
#'  and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the basic
#' triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle
#' (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:arc-density-PE;textual}{pcds}).
#' Hence, standard basic triangle is useful for simulation studies
#' under the uniformity hypothesis.
#'
#' @param n A positive integer representing the number of uniform points
#' to be generated in the standard basic triangle.
#' @param c1,c2 Positive real numbers representing the top vertex
#' in standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and
#' \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{tess.points}{The vertices of the support
#' of the uniformly generated points,
#' it is the standard basic triangle \eqn{T_b} for this function}
#' \item{gen.points}{The output set of generated points uniformly
#' in the standard basic triangle}
#' \item{out.region}{The outer region which contains the support region,
#' \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern
#' from which points are to be generated}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number
#' of vertices of the support points (here it is 3).}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the support, Tb}
#'
#' @seealso \code{\link{runif.std.tri}}, \code{\link{runif.tri}},
#' and \code{\link{runif.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C);
#' n<-100
#'
#' set.seed(1)
#' runif.basic.tri(1,c1,c2)
#' Xdt<-runif.basic.tri(n,c1,c2)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-runif.basic.tri(n,c1,c2)$g
#'
#' Xlim<-range(Tb[,1])
#' Ylim<-range(Tb[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,xlab="",ylab="",xlim=Xlim+xd*c(-.01,.01),
#' ylim=Ylim+yd*c(-.01,.01),type="n")
#' polygon(Tb)
#' points(Xp)
#' }
#'
#' @export runif.basic.tri
runif.basic.tri <- function(n,c1,c2)
{
  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  X <-matrix(0,n,2);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,0,1); y<-runif(1,0,c2);
    if (y < (c2/c1)*x && y < c2*(x-1)/(c1-1))
    {X[i,] <- c(x,y);
    ct<-1;
    }
  }
  }
  A<-c(0,0); B<-c(1,0); C<-c(c1,c2)
  Tb<-rbind(A,B,C)  #standard basic triangle
  row.names(Tb) = NULL

  Cvec= paste(round(c(c1,c2),2), collapse=",");
  typ<-paste("Uniform Distribution in the Standard Basic Triangle with vertices (0,0), (1,0), and (",Cvec,")",sep="")
  main.txt<-paste("Uniform Points in the Standard Basic Triangle\n with vertices (0,0), (1,0), and (",Cvec,")",sep="")
  txt<-paste(n, " uniform points in the standard basic triangle with vertices (0,0), (1,0), and (",Cvec,")",sep="")

  npts<-c(n,3)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], "is the number of uniform points \n",names(npts)[2], "is the number of vertices of the support region")

  res<-list(
    type=typ,
    gen.points=X,
    #uniformly generated points in the standard basic triangle, Tb
    tess.points=Tb,
    #tessellation points whose convex hull constitutes the support of the uniform points
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    out.region=NULL, #outer region for Tb
    num.points=npts,  txt4pnts=txt4pnts,
    xlimit=range(Tb[,1]),
    ylimit=range(Tb[,2])
  )

  class(res)<-"Uniform"
  res$call <-match.call()
  res
} #end of the function
#'


#################################################################

#' @title Generation of Uniform Points in a Triangle
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly in a given triangle, \code{tri}
#'
#' @param n A positive integer representing the number of uniform points
#' to be generated in the triangle.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{tess.points}{The vertices of the support of
#' the uniformly generated points, it is the triangle
#' \code{tri} for this function}
#' \item{gen.points}{The output set of generated points uniformly
#' in the triangle, \code{tri}.}
#' \item{out.region}{The outer region which contains the support region,
#' \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern
#' from which points are to be generated}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number
#' of vertices of the support points (here it is 3).}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the support, \code{tri}}
#'
#' @seealso \code{\link{runif.std.tri}}, \code{\link{runif.basic.tri}},
#' and \code{\link{runif.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-100
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C)
#'
#' Xdt<-runif.tri(n,Tr)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-Xdt$g
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' plot(Tr,pch=".",xlab="",ylab="",main="Uniform Points in One Triangle",
#'      xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp)
#' xc<-Tr[,1]+c(-.02,.02,.02)
#' yc<-Tr[,2]+c(.02,.02,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#' }
#'
#' @export runif.tri
runif.tri <- function(n,tri)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  x1<-min(tri[,1]); x2<-max(tri[,1])
  y1<-min(tri[,2]); y2<-max(tri[,2])

  X <-matrix(0,n,2);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,x1,x2); y<-runif(1,y1,y2);
    if (in.triangle(c(x,y),tri)$in.tri==TRUE)
    {X[i,]<-c(x,y);
    ct<-1;
    }
  }
  }

  Avec= paste(round(tri[1,],2), collapse=",");
  Bvec= paste(round(tri[2,],2), collapse=",");
  Cvec= paste(round(tri[3,],2), collapse=",");
  main.txt<-paste("Uniform Points in the Triangle \n with Vertices (",Avec,"), (",Bvec,") and (",Cvec,")",sep="")
  txt<-paste(n, " uniform points in the triangle with vertices (",Avec,"), (",Bvec,") and (",Cvec,")",sep="")
  typ<-paste("Uniform Distribution in the Triangle with Vertices (",Avec,"), (",Bvec,"), and (",Cvec,")",sep="")

  npts<-c(n,3)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], " : the number of uniform points\n",names(npts)[2], " : the number of vertices of the support region", sep="")

  res<-list(
    type=typ,
    mtitle=main.txt,
    gen.points=X, #uniformly generated points in the triangle, tri
    tess.points=tri,
    #tessellation points whose convex hull constitutes the support of the uniform points
    desc.pat=txt, #description of the pattern
    out.region=NULL, #outer region for tri
    num.points=npts,  txt4pnts=txt4pnts,
    xlimit=range(tri[,1]),
    ylimit=range(tri[,2])
  )

  class(res)<-"Uniform"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points segregated (in a Type I fashion)
#' from the vertices of \eqn{T_e}
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly
#' in the standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' under the type I segregation alternative for \code{eps}
#' in \eqn{(0,\sqrt{3}/3=0.5773503]}.
#'
#' In the type I segregation, the triangular forbidden regions
#' around the vertices are determined by
#' the parameter \code{eps}
#' which serves as the height of these triangles
#' (see examples for a sample plot.)
#'
#' See also (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param eps A positive real number
#' representing the parameter of type I segregation (which is the
#' height of the triangular forbidden regions around the vertices).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{The exclusion parameter, \code{eps},
#' of the segregation pattern, which is the height
#' of the triangular forbidden regions around the vertices }
#' \item{ref.points}{The input set of points \code{Y};
#' reference points, i.e., points
#' from which generated points are segregated
#' (i.e., vertices of \eqn{T_e}).}
#' \item{gen.points}{The output set of generated points segregated
#' from \code{Y} points (i.e., vertices of \eqn{T_e}).}
#' \item{tri.Y}{Logical output for triangulation
#' based on \code{Y} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Y} points is
#' to be implemented (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and
#' the number of reference (i.e., \code{Y}) points, which is 3 here.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points,
#' which are the vertices of \eqn{T_e} here.}
#'
#' @seealso \code{\link{rseg.circular}}, \code{\link{rassoc.circular}},
#' \code{\link{rsegII.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' n<-100
#' eps<-.3  #try also .15, .5, .75
#'
#' set.seed(1)
#' Xdt<-rseg.std.tri(n,eps)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' Xp<-Xdt$gen.points
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",
#' main="Type I segregation in the \n standard equilateral triangle",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp)
#'
#' #The support for the Type I segregation alternative
#' sr<-eps/(sqrt(3)/2)
#' C1<-C+sr*(A-C); C2<-C+sr*(B-C)
#' A1<-A+sr*(B-A); A2<-A+sr*(C-A)
#' B1<-B+sr*(A-B); B2<-B+sr*(C-B)
#' supp<-rbind(A1,B1,B2,C2,C1,A2)
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",
#' main="Support of the Type I Segregation",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' if (sr<=.5)
#' {
#'   polygon(Te)
#'   polygon(supp,col=5)
#' } else
#' {
#'   polygon(Te,col=5,lwd=2.5)
#'   polygon(rbind(A,A1,A2),col=0,border=NA)
#'   polygon(rbind(B,B1,B2),col=0,border=NA)
#'   polygon(rbind(C,C1,C2),col=0,border=NA)
#' }
#' points(Xp)
#' }
#'
#' @export rseg.std.tri
rseg.std.tri <- function(n,eps)
{
  if (!is.point(eps,1) || eps<=0 || eps>=sqrt(3)/3)
  {stop('eps must be a scalar in (0,sqrt(3)/3=0.5773503)')}

  Y<-rbind(c(0,0),c(1,0),c(1/2,sqrt(3)/2))  #std eq triangle
  ny<-nrow(Y)

  X <-matrix(0,n,2);
  if (eps < 0.4330127 )  #(eps < sqrt(3)/4 )
  {
    for (i in 1:n)
    {ct<-0;
    while (ct==0)
    {
      x<-runif(1,0.5773503*eps,1-0.5773503*eps); y<-runif(1,0,0.8660254-eps);
      if (y < 1.732050808*x && y < 1.732050808-1.732050808*x &&
          y > 1.732051*x+2*eps-1.732051 && y > -1.732051*x+2*eps)
      {X[i,]<-c(x,y);
      ct<-1;
      }
    }
    }
  }
  else
  {
    for (i in 1:n)
    {ct<-0;
    while (ct==0)
    {
      x<-runif(1,-.5+1.732050809*eps,1.5-1.732050809*eps); y<-runif(1,-.866025404+2*eps,.866025404-eps);
      if (y > 2*eps-1.732050808*x && y > 2*eps-1.732050808+1.732050808*x)
      {X[i,]<-c(x,y);
      ct<-1;
      }
    }
    }
  }

  param<-eps
  names(param)<-"exclusion parameter"
  typ<-paste("Type I Segregation of ",n, " points from vertices of the standard equilateral triangle with exclusion parameter eps = ",eps,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Type I Segregation of One Class from Vertices of the Standard Equilateral Triangle"
  main.txt<-paste("Type I Segregation in the Standard Equilateral Triangle \n with Exclusion Parameter eps = ",eps,sep="")

  res<-list(
    type=typ,
    parameters=param,
    gen.points=X, #generated points segregated from Y points (vertices of std eq triangle)
    ref.points=Y, #reference points, i.e., points from which generated points are segregated (vertices of std eq triangle)
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    tri.Y=TRUE,
    num.points=npts,
    xlimit=range(Y[,1]),
    ylimit=range(Y[,2])
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points segregated (in a Type I fashion)
#' from the vertices of a triangle
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly in the support
#' for Type I segregation in a given triangle, \code{tri}.
#'
#' \code{delta} is the parameter of segregation (that is,
#' \eqn{\delta 100} \% of the area around each vertex
#' in the triangle is forbidden for point generation).
#' \code{delta} corresponds to \code{eps} in the
#' standard equilateral triangle \eqn{T_e} as \eqn{delta=4 eps^2/3}
#' (see \code{rseg.std.tri}  function).
#'
#' See (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the segregation pattern.
#'
#' @param n A positive integer representing the number of points
#' to be generated from the segregation pattern
#' in the triangle, \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param delta A positive real number in \eqn{(0,4/9)}.
#' \code{delta} is the parameter of segregation (that is,
#' \eqn{\delta 100} \% area around each vertex
#' in each Delaunay triangle is forbidden for point generation).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title
#' for the plot of the point pattern}
#' \item{parameters}{Exclusion parameter, \code{delta},
#' of the Type I segregation pattern.
#' \code{delta} is in \eqn{(0,4/9)}
#' \eqn{\delta 100} \% area around each vertex
#' in the triangle \code{tri} is forbidden for point generation.}
#' \item{ref.points}{The input set of points, i.e., vertices of \code{tri};
#' reference points, i.e., points
#' from which generated points are segregated.}
#' \item{gen.points}{The output set of generated points segregated
#' from the vertices of \code{tri}.}
#' \item{tri.Y}{Logical output,
#' if \code{TRUE} the triangle \code{tri} is also plotted when the
#' corresponding plot function from the \code{Patterns} object is called.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points
#' and the number of reference (i.e., vertex of \code{tri},
#' which is 3 here).}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points,
#' which are the vertices of the triangle \code{tri}}
#'
#' @seealso \code{\link{rassoc.tri}}, \code{\link{rseg.std.tri}},
#' \code{\link{rsegII.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-100
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C)
#' del<-.4
#'
#' Xdt<-rseg.tri(n,Tr,del)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-Xdt$g
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",
#' main="Points from Type I Segregation \n in one Triangle",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp)
#' xc<-Tr[,1]+c(-.02,.02,.02)
#' yc<-Tr[,2]+c(.02,.02,.03)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
rseg.tri <- function(n,tri,delta)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  if (!is.point(delta,1) || delta<=0 || delta>=1)
  {stop('delta must be a scalar in (0,1)')}

  cnt<-1;ind<-0
  Xp<-NULL
  while (cnt <= n)
  {
    pnt<-as.vector(runif.tri(1,tri)$g)
    seg.tri<-seg.tri.support(delta,tri)
    if (in.triangle(pnt,seg.tri,boundary=TRUE)$in.tri)
    {
      Xp<-rbind(Xp,pnt);
      cnt<-cnt+1
    }
  }

  param<-delta
  names(param)<-"exclusion parameter"

  Avec= paste(round(tri[1,],2), collapse=",");
  Bvec= paste(round(tri[2,],2), collapse=",");
  Cvec= paste(round(tri[3,],2), collapse=",");
  typ<-paste("Type I Segregation of ",n, " points in the triangle with vertices (",Avec,"), ",
             "(",Bvec,")"," and ","(",Cvec,") and exclusion parameter delta = ",delta,sep="")

  npts<-c(n,3)
  names(npts)<-c("nx","ny")

  txt<-"Type I Segregation of One Class from Vertices of the Triangle"
  main.txt<-paste("Type I Segregation in the One-Triangle Case \n with Exclusion Parameter delta = ",delta,sep="")

  res<-list(
    type=typ,
    parameters=param,
    gen.points=Xp,
    #generated points segregated from Y points (vertices of the triangle, tri)
    ref.points=tri,
    #reference points, i.e., points from which generated points are segregated (vertices of the triangle)
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    tri.Y=TRUE,
    num.points=npts,
    xlimit=range(tri[,1]),
    ylimit=range(tri[,2])
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points segregated (in a Type II fashion)
#' from the vertices of \eqn{T_e}
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly
#' in the standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' under the type II segregation alternative
#' for \code{eps} in \eqn{(0,\sqrt{3}/6=0.2886751]}.
#'
#' In the type II segregation, the annular forbidden regions
#' around the edges are determined by
#' the parameter \code{eps} which is the distance from the interior triangle
#' (i.e., support for the segregation)
#' to \eqn{T_e} (see examples for a sample plot.)
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param eps A positive real number
#' representing the parameter of type II segregation (which is the
#' distance from the interior triangle points to the boundary of \eqn{T_e}).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{The exclusion parameter,
#' \code{eps}, of the segregation pattern,
#' which is the distance from the interior triangle to \eqn{T_e}}
#' \item{ref.points}{The input set of points \code{Y};
#' reference points, i.e., points
#' from which generated points are segregated (i.e., vertices of \eqn{T_e}).}
#' \item{gen.points}{The output set of generated points
#' segregated from \code{Y} points
#' (i.e., vertices of \eqn{T_e}).}
#' \item{tri.Y}{Logical output for triangulation
#' based on \code{Y} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Y} points is to be implemented
#' (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and
#' the number of reference (i.e., \code{Y}) points, which is 3 here.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points,
#' which are the vertices of \eqn{T_e} here}
#'
#' @seealso \code{\link{rseg.circular}}, \code{\link{rassoc.circular}},
#' \code{\link{rseg.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' n<-10  #try also n<-20 or n<-100 or 1000
#' eps<-.15  #try also .2
#'
#' set.seed(1)
#' Xdt<-rsegII.std.tri(n,eps)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' Xp<-Xdt$gen.points
#'
#' plot(Te,pch=".",xlab="",ylab="",
#' main="Type II segregation in the \n standard equilateral triangle",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp)
#'
#' #The support for the Type II segregation alternative
#' C1<-c(1/2,sqrt(3)/2-2*eps);
#' A1<-c(eps*sqrt(3),eps); B1<-c(1-eps*sqrt(3),eps);
#' supp<-rbind(A1,B1,C1)
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",
#' main="Support of the Type II Segregation",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#'   polygon(Te)
#'   polygon(supp,col=5)
#' points(Xp)
#'
#' @export rsegII.std.tri
rsegII.std.tri <- function(n,eps)
{
  if (!is.point(eps,1) || eps<=0 || eps>=sqrt(3)/6)
  {stop('eps must be a scalar in (0,sqrt(3)/6=0.2886751)')}

  Y<-rbind(c(0,0),c(1,0),c(1/2,sqrt(3)/2))  #std eq triangle
  ny<-nrow(Y)

  X <-matrix(0,n,2);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,1.732050808*eps,1-1.732050808*eps);
    y<-runif(1,eps,0.8660254-2*eps);
    if (y < -2*eps+1.732050808*x &&
        y < -2*eps+1.732050808-1.732050808*x)
    {X[i,]<-c(x,y);
    ct<-1;
    }
  }
  }

  param<-eps
  names(param)<-"exclusion parameter"
  typ<-paste("Type II Segregation of ",n, " points from edges of the standard equilateral triangle with exclusion parameter eps = ",eps,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Type II Segregation of One Class from Edges of the Standard Equilateral Triangle"
  main.txt<-paste("Type II Segregation in the Standard Equilateral Triangle \n with Exclusion Parameter eps = ",eps,sep="")

  res<-list(
    type=typ,
    parameters=param,
    gen.points=X,
    #generated points segregated from Y points (vertices of std eq triangle)
    ref.points=Y,
    #reference points, i.e., points from which generated points are segregated (vertices of std eq triangle)
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    tri.Y=TRUE,
    num.points=npts,
    xlimit=range(Y[,1]),
    ylimit=range(Y[,2])
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points associated (in a Type I fashion)
#' with the vertices of \eqn{T_e}
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly
#' in the standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' under the type I association alternative
#' for \code{eps} in \eqn{(0,\sqrt{3}/3=0.5773503]}.
#' The allowed triangular regions around the vertices are determined
#' by the parameter \code{eps}.
#'
#' In the type I association, the triangular support regions
#' around the vertices are determined by
#' the parameter \code{eps}
#' where \eqn{\sqrt{3}/3}-\code{eps} serves as the height of these triangles
#' (see examples for a sample plot.)
#'
#' See also (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param eps A positive real number
#' representing the parameter of type I association
#' (where \eqn{\sqrt{3}/3}-\code{eps}
#' serves as the height of the triangular support regions
#' around the vertices).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{The attraction parameter of the association pattern,
#' \code{eps},
#' where \eqn{\sqrt{3}/3}-\code{eps} serves
#' as the height of the triangular support regions around the vertices}
#' \item{ref.points}{The input set of points \code{Y};
#' reference points, i.e., points
#' with which generated points are associated (i.e., vertices of \eqn{T_e}).}
#' \item{gen.points}{The output set of generated points
#' associated with \code{Y} points (i.e., vertices of \eqn{T_e}).}
#' \item{tri.Y}{Logical output for triangulation
#' based on \code{Y} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Y} points is
#' to be implemented (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern.}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points
#' and the number of reference (i.e., \code{Y}) points.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points, which are the
#' vertices of \eqn{T_e} here}
#'
#' @seealso \code{\link{rseg.circular}}, \code{\link{rassoc.circular}},
#' \code{\link{rsegII.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' n<-100  #try also n<-20 or n<-100 or 1000
#' eps<-.25  #try also .15, .5, .75
#'
#' set.seed(1)
#' Xdt<-rassoc.std.tri(n,eps)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' Xp<-Xdt$gen.points
#' plot(Te,pch=".",xlab="",ylab="",
#' main="Type I association in the \n standard equilateral triangle",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp)
#'
#' #The support for the Type I association alternative
#' sr<-(sqrt(3)/3-eps)/(sqrt(3)/2)
#' C1<-C+sr*(A-C); C2<-C+sr*(B-C)
#' A1<-A+sr*(B-A); A2<-A+sr*(C-A)
#' B1<-B+sr*(A-B); B2<-B+sr*(C-B)
#' supp<-rbind(A1,B1,B2,C2,C1,A2)
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",
#' main="Support of the Type I Association",
#' xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' if (sr<=.5)
#' {
#'   polygon(Te,col=5)
#'   polygon(supp,col=0)
#' } else
#' {
#'   polygon(Te,col=0,lwd=2.5)
#'   polygon(rbind(A,A1,A2),col=5,border=NA)
#'   polygon(rbind(B,B1,B2),col=5,border=NA)
#'   polygon(rbind(C,C1,C2),col=5,border=NA)
#' }
#' points(Xp)
#' }
#'
#' @export rassoc.std.tri
rassoc.std.tri <- function(n,eps)
{
  if (!is.point(eps,1) || eps<=0 || eps>=sqrt(3)/3 )
  {stop('eps must be a scalar in (0,sqrt(3)/3=0.5773503)')}

  Y<-rbind(c(0,0),c(1,0),c(1/2,sqrt(3)/2))  #std eq triangle
  ny<-nrow(Y)

  X<-matrix(0,n,2);
  Eps<-0.5773503-eps; #sqrt(3)/3-eps
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,0,1); y<-runif(1,0,0.8660254);
    if ( (y < 1.732050808*x && y < 1.732050808-1.732050808*x) &&
         !(y > 1.732051*x+2*Eps-1.732051 &
           y > -1.732051*x+2*Eps & y < 0.8660254-Eps) )
    {X[i,]<-c(x,y);
    ct<-1;
    }
  }
  }

  param<-eps
  names(param)<-"attraction parameter"
  typ<-paste("Type I Association of ",n, " points with vertices of the standard equilateral triangle with attraction parameter eps = ",eps,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Type I Association of One Class with Vertices of the Standard Equilateral Triangle"
  main.txt<-paste("Type I Association in the Standard Equilateral Triangle \n with Attraction Parameter eps =",eps,sep="")

  res<-list(
    type=typ,
    parameters=param,
    gen.points=X,
    #generated points associated with Y points (vertices of std eq triangle)
    ref.points=Y,
    #reference points, i.e., points with which generated points are associated (vertices of std eq triangle)
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    tri.Y=TRUE,
    num.points=npts,
    xlimit=range(Y[,1]),
    ylimit=range(Y[,2])
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points associated (in a Type I fashion)
#' with the vertices of a triangle
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly in the support
#'  for Type I association in a given triangle, \code{tri}.
#' \code{delta} is the parameter of association
#' (that is, only \eqn{\delta 100} \% area around each vertex in
#' the triangle is allowed for point generation).
#' \code{delta} corresponds to \code{eps}
#' in the standard equilateral triangle
#' \eqn{T_e} as \eqn{delta=4eps^2/3}
#' (see \code{rseg.std.tri} function).
#'
#' See (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}) for more on
#' the association pattern.
#'
#' @param n A positive integer representing the number of points
#' to be generated from the association pattern
#' in the triangle, \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param delta A positive real number in \eqn{(0,4/9)}.
#' \code{delta} is the parameter of association
#' (that is, only \eqn{\delta 100} \% area around each vertex
#' in the triangle is allowed for point generation).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title
#' for the plot of the point pattern}
#' \item{parameters}{Attraction parameter, \code{delta},
#' of the Type I association pattern.
#' \code{delta} is in \eqn{(0,4/9)}
#' only \eqn{\delta 100} \% of the area around each vertex
#' in the triangle \code{tri}
#' is allowed for point generation.}
#' \item{ref.points}{The input set of points,
#' i.e., vertices of \code{tri};
#' reference points, i.e., points
#' with which generated points are associated.}
#' \item{gen.points}{The output set of generated points
#' associated with the vertices of \code{tri}.}
#' \item{tri.Y}{Logical output, \code{TRUE} if triangulation
#' based on \code{Yp} points should be implemented.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points
#' and the number of reference (i.e., \code{Yp}) points.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}- and
#' \eqn{y}-coordinates of the reference points, which are the
#' \code{Yp} points}
#'
#' @seealso \code{\link{rseg.tri}}, \code{\link{rassoc.std.tri}},
#' \code{\link{rassocII.std.tri}}, and \code{\link{rassoc.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-100
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C)
#' del<-.4
#'
#' Xdt<-rassoc.tri(n,Tr,del)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-Xdt$g
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",
#' main="Points from Type I Association \n in one Triangle",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp)
#' xc<-Tr[,1]+c(-.02,.02,.02)
#' yc<-Tr[,2]+c(.02,.02,.03)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
rassoc.tri <- function(n,tri,delta)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  if (!is.point(delta,1) || delta<=0 || delta>=1)
  {stop('delta must be a scalar in (0,1)')}

  cnt<-1;ind<-0
  Xp<-NULL
  while (cnt <= n)
  {
    pnt<-as.vector(runif.tri(1,tri)$g)
    seg.tri<-seg.tri.support(delta,tri)
    if (!in.triangle(pnt,seg.tri)$in.tri)
    {
      Xp<-rbind(Xp,pnt);
      cnt<-cnt+1
    }
  }

  param<-delta
  names(param)<-"attraction parameter"

  Avec= paste(round(tri[1,],2), collapse=",");
  Bvec= paste(round(tri[2,],2), collapse=",");
  Cvec= paste(round(tri[3,],2), collapse=",");
  typ<-paste("Type I Association of ",n, " points in the triangle with vertices (",Avec,"), ",
             "(",Bvec,")"," and ","(",Cvec,") with attraction parameter delta = ",delta,sep="")

  npts<-c(n,3)
  names(npts)<-c("nx","ny")

  txt<-"Type I Association of One Class with Vertices of the Triangle"
  main.txt<-paste("Type I Association in the One-Triangle Case \n with Attraction Parameter delta = ",delta,sep="")

  res<-list(
    type=typ,
    parameters=param,
    gen.points=Xp, #generated points associated with Y points (vertices of the triangle, tri)
    ref.points=tri, #reference points, i.e., points with which generated points are associated (vertices of the triangle)
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    tri.Y=TRUE,
    num.points=npts,
    xlimit=range(tri[,1]),
    ylimit=range(tri[,2])
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of points associated (in a Type II fashion)
#' with the edges of \eqn{T_e}
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly
#' in the standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' under the type II association alternative for \code{eps}
#' in \eqn{(0,\sqrt{3}/6=0.2886751]}.
#'
#' In the type II association, the annular allowed regions
#' around the edges are determined by
#' the parameter \code{eps}
#' where \eqn{\sqrt{3}/6}-\code{eps} is the distance
#' from the interior triangle
#' (i.e., forbidden region for association) to \eqn{T_e}
#' (see examples for a sample plot.)
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param eps A positive real number
#' representing the parameter of type II association
#' (where \eqn{\sqrt{3}/6}-\code{eps}
#' is the distance from the interior triangle distance
#' from the interior triangle to \eqn{T_e}).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the point pattern}
#' \item{mtitle}{The \code{"main"} title
#' for the plot of the point pattern}
#' \item{parameters}{The attraction parameter, \code{eps},
#' of the association pattern,
#' where \eqn{\sqrt{3}/6}-\code{eps}
#' is the distance from the interior triangle to \eqn{T_e}}
#' \item{ref.points}{The input set of points \code{Y};
#' reference points, i.e.,
#' points with which generated points are associated
#' (i.e., vertices of \eqn{T_e}).}
#' \item{gen.points}{The output set of generated points associated
#' with \code{Y} points (i.e., edges of \eqn{T_e}).}
#' \item{tri.Y}{Logical output for triangulation
#' based on \code{Y} points should be implemented or not.
#' if \code{TRUE} triangulation based on \code{Y} points is
#' to be implemented (default is set to \code{FALSE}).}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and
#' the number of reference (i.e., \code{Y}) points,
#' which is 3 here.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points,
#' which are the vertices of \eqn{T_e} here.}
#'
#' @seealso \code{\link{rseg.circular}}, \code{\link{rassoc.circular}},
#' \code{\link{rsegII.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' n<-100  #try also n<-20 or n<-100 or 1000
#' eps<-.2  #try also .25, .1
#'
#' set.seed(1)
#' Xdt<-rassocII.std.tri(n,eps)
#' Xdt
#' summary(Xdt)
#' plot(Xdt,asp=1)
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' Xp<-Xdt$gen.points
#' plot(Te,pch=".",xlab="",ylab="",
#' main="Type II association in the \n standard equilateral triangle",
#'      xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp)
#'
#' #The support for the Type II association alternative
#' A1<-c(1/2-eps*sqrt(3),sqrt(3)/6-eps);
#' B1<-c(1/2+eps*sqrt(3),sqrt(3)/6-eps);
#' C1<-c(1/2,sqrt(3)/6+2*eps);
#' supp<-rbind(A1,B1,C1)
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",
#' main="Support of the Type II Association",
#' xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te,col=5)
#' polygon(supp,col=0)
#' points(Xp)
#' }
#'
#' @export rassocII.std.tri
rassocII.std.tri <- function(n,eps)
{
  if (!is.point(eps,1) || eps<=0 || eps>=sqrt(3)/6 )
  {stop('eps must be a scalar in (0,sqrt(3)/6=0.2886751)')}

  Y<-rbind(c(0,0),c(1,0),c(1/2,sqrt(3)/2))  #std eq triangle
  ny<-nrow(Y)

  X <-matrix(0,n,2);
  Eps<-0.2886751-eps; #sqrt(3)/6-eps;
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,0,1); y<-runif(1,0,0.8660254);
    if ( (y < 1.732050808*x && y < 1.732050808-1.732050808*x) &&
         !(y > Eps & y < -2*Eps+1.732050808*x  &
           y < -2*Eps+1.732050808-1.732050808*x) )
    {X[i,]<-c(x,y);
    ct<-1;
    }
  }
  }

  param<-eps
  names(param)<-"attraction parameter"
  typ<-paste("Type II Association of ",n, " points with edges of the standard equilateral triangle with attraction parameter eps = ",eps,sep="")

  npts<-c(n,ny)
  names(npts)<-c("nx","ny")

  txt<-"Type II Association of One Class with Edges of the Standard Equilateral Triangle"
  main.txt<-paste("Type II Association Pattern in the Standard Equilateral Triangle \n with Attraction Parameter eps = ",eps,sep="")
  res<-list(
    type=typ,
    parameters=param,
    gen.points=X,
    #generated points associated with Y points (edges of std eq triangle)
    ref.points=Y,
    #reference points, i.e., points with which generated points are associated (edges eq triangle)
    desc.pat=txt, #description of the pattern
    mtitle=main.txt,
    tri.Y=TRUE,
    num.points=npts,
    xlimit=range(Y[,1]),
    ylimit=range(Y[,2])
  )

  class(res)<-"Patterns"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#functions generated for general triangles

#' @title Generation of Uniform Points in the Convex Hull of Points
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly
#' in the Convex Hull of set of points, \code{Yp}.
#' That is, generates uniformly in each of the triangles
#' in the Delaunay triangulation of \code{Yp}, i.e.,
#' in the multiple triangles partitioning the convex hull of \code{Yp}.
#'
#' If \code{Yp} consists only of 3 points,
#' then the function behaves like the
#' function \code{\link{runif.tri}}.
#'
#' \code{DTmesh} is the Delaunay triangulation of \code{Yp},
#' default is \code{DTmesh=NULL}.
#' \code{DTmesh} yields triangulation nodes with neighbours
#' (result of \code{\link[interp]{tri.mesh}} function
#' from \code{interp} package).
#'
#' See (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param n A positive integer
#' representing the number of uniform points to be generated
#' in the convex hull of the point set \code{Yp}.
#' @param Yp A set of 2D points
#' whose convex hull is the support of the uniform points to be generated.
#' @param DTmesh Triangulation nodes with neighbours
#' (result of \code{\link[interp]{tri.mesh}} function from
#' \code{interp} package).
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title
#' for the plot of the point pattern}
#' \item{tess.points}{The points
#' which constitute the vertices of the triangulation and
#' whose convex hull determines the support of the generated points.}
#' \item{gen.points}{The output set of generated points uniformly
#' in the convex hull of \code{Yp}}
#' \item{out.region}{The outer region which contains the support region,
#' \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern
#' from which points are to be generated}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and
#' the number of vertices in the triangulation
#' (i.e., size of \code{Yp}) points.}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the points in Yp}
#'
#' @seealso \code{\link{runif.tri}}, \code{\link{runif.std.tri}},
#' and \code{\link{runif.basic.tri}},
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-4;  #try also nx<-1000; ny<-10;
#' set.seed(1)
#' Yp<-cbind(runif(ny,0,10),runif(ny,0,10))
#'
#' Xdt<-runif.multi.tri(nx,Yp)
#' #data under CSR in the convex hull of Ypoints
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-Xdt$g
#' #or use
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
#' #Delaunay triangulation based on Y points
#' Xp<-runif.multi.tri(nx,Yp,DTY)$g
#' #data under CSR in the convex hull of Ypoints
#'
#' Xlim<-range(Yp[,1])
#' Ylim<-range(Yp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' #plot of the data in the convex hull of Y points together with the Delaunay triangulation
#' plot(Xp, xlab=" ", ylab=" ",
#' main="Uniform Points in Convex Hull of Y Points",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),type="n")
#' interp::plot.triSht(DTY, add=TRUE,
#' do.points = TRUE,pch=16,col="blue")
#' points(Xp,pch=".",cex=3)
#'
#' Yp<-rbind(c(.3,.2),c(.4,.5),c(.14,.15))
#' runif.multi.tri(nx,Yp)
#' }
#'
#' @export runif.multi.tri
runif.multi.tri <- function(n,Yp,DTmesh=NULL)
{
  Yp<-as.matrix(Yp)
  if (!is.numeric(Yp) || ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be numeric and of dimension kx2 with k>=3')}

  ny<-nrow(Yp)
  if (ny==3)
  {
    res<-runif.tri(n,Yp)
  } else
  {
  if (is.null(DTmesh))
  {DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")}
    #Delaunay triangulation

  cnt<-1
  Xp<-NULL
  while (cnt <= n)
  {
    lx<-min(Yp[,1]); rx<-max(Yp[,1])
    ly<-min(Yp[,2]); ry<-max(Yp[,2])
    x1<-runif(1,lx,rx); y1<-runif(1,ly,ry)
    if (interp::in.convex.hull(DTmesh,x1,y1,strict=FALSE))
    {
      Xp<-rbind(Xp,c(x1,y1));
      cnt<-cnt+1
    }
  }
  typ<-"Uniform distribution in the convex hull of Yp points"
  main.txt<-"Uniform Points in Convex Hull of Yp Points"
  txt<-paste(n, " uniform points in the convex hull of Yp points",sep="")

  ny<-nrow(Yp)
  npts<-c(n,ny)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], " : the number of Uniform points\n",
                  names(npts)[2], " : the number of points whose convex hull determines the support",sep="")

  res<-list(
    type=typ,
    mtitle=main.txt,
    gen.points=Xp,
    #uniformly generated points in the triangle, tri
    tess.points=Yp,
    #The points that determine the support of the uniform points
    desc.pat=txt, #description of the pattern
    out.region=NULL, #outer region for convex hull
    num.points=npts,  txt4pnts=txt4pnts,
    xlimit=range(Yp[,1]),
    ylimit=range(Yp[,2])
  )

  class(res)<-"Uniform"
  res$call <-match.call()
  }

  res
} #end of the function
#'

#################################################################

#' @title The auxiliary triangle to define the support of type I segregation
#'
#' @description Returns the triangle whose intersection
#' with a general triangle gives the support for
#' type I segregation given the \code{delta}
#' (i.e., \eqn{\delta 100} \% area of a triangle around the
#' vertices is chopped off).
#' See the plot in the examples.
#'
#' Caveat: the vertices of this triangle may be
#' outside the triangle, \code{tri}, depending on the value of
#' \code{delta} (i.e., for small values of \code{delta}).
#'
#' @param delta A positive real number between 0 and 1
#' that determines the percentage of area of the triangle
#' around the vertices forbidden for point generation.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#'
#' @return the vertices of the triangle (stacked row-wise)
#' whose intersection with a general triangle
#' gives the support for type I segregation for the given \code{delta}
#'
#' @seealso \code{\link{rseg.std.tri}} and \code{\link{rseg.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #for a general triangle
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' delta<-.3  #try also .5,.75,.85
#' Tseg<-seg.tri.support(delta,Tr)
#'
#' Xlim<-range(Tr[,1],Tseg[,1])
#' Ylim<-range(Tr[,2],Tseg[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' par(pty="s")
#' plot(Tr,pch=".",xlab="",ylab="",
#' main="segregation support is the intersection\n of these two triangles",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' polygon(Tseg,lty=2)
#'
#' txt<-rbind(Tr,Tseg)
#' xc<-txt[,1]+c(-.03,.03,.03,.06,.04,-.04)
#' yc<-txt[,2]+c(.02,.02,.04,-.03,0,0)
#' txt.str<-c("A","B","C","T1","T2","T3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
seg.tri.support <- function(delta,tri)
{
  if (!is.point(delta,1) || delta<=0 || delta>=1)
  {stop('delta must be a scalar in (0,1)')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate')}

  if (delta <= 3/4)
  {k<-sqrt(delta/3);}
  else
  {k<-(2-sqrt(1-delta))/3}
  A<-tri[1,]; a1<-A[1]; a2<-A[2];
  B<-tri[2,]; b1<-B[1]; b2<-B[2];
  C<-tri[3,]; c1<-C[1]; c2<-C[2];

  T1<-c(2*k*c1-c1+a1-b1*k+b1-a1*k, -b2*k-a2*k+a2+2*c2*k-c2+b2);
  T2<-c(-k*c1+c1-a1+2*a1*k-b1*k+b1, 2*a2*k-b2*k+b2-c2*k+c2-a2);
  T3<-c(-k*c1+c1-b1-a1*k+a1+2*b1*k, -a2*k+2*b2*k+a2-c2*k+c2-b2);

  TRI<-rbind(T1,T2,T3)
  row.names(TRI)<-c()
  TRI
} #end of the function
#'

#################################################################

#' @title Generation of points segregated (in a Type I fashion)
#' from a given set of points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly in the support
#' for Type I segregation in the convex hull of
#' set of points, \code{Yp}.
#'
#' \code{delta} is the parameter of segregation
#' (that is, \eqn{\delta 100} \% of the area around each vertex
#' in each Delaunay
#' triangle is forbidden for point generation).
#' \code{delta} corresponds to \code{eps}
#' in the standard equilateral triangle
#' \eqn{T_e} as \eqn{delta=4 eps ^2/3} (see \code{rseg.std.tri} function).
#'
#' If \code{Yp} consists only of 3 points,
#' then the function behaves like the
#' function \code{\link{rseg.tri}}.
#'
#' \code{DTmesh} must be the Delaunay triangulation of \code{Yp}
#' and \code{DTr} must be the corresponding Delaunay triangles
#' (both \code{DTmesh} and \code{DTr} are \code{NULL} by default).
#' If \code{NULL}, \code{DTmesh} is computed
#' via \code{\link[interp]{tri.mesh}}
#' and \code{DTr} is computed via \code{\link[interp]{triangles}}
#' function in \code{interp} package.
#'
#' \code{\link[interp]{tri.mesh}} function yields the triangulation nodes
#' with their neighbours,
#' and creates a triangulation object,
#' and \code{\link[interp]{triangles}} function
#' yields a triangulation data structure
#' from the triangulation object created
#' by \code{\link[interp]{tri.mesh}} (the first three columns are
#' the vertex indices of the Delaunay triangles.)
#'
#' See (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the segregation pattern.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param Yp A set of 2D points
#' from which Delaunay triangulation is constructed.
#' @param delta A positive real number in \eqn{(0,4/9)}.
#' \code{delta} is the parameter of segregation
#' (that is, \eqn{\delta 100} %
#' area around each vertex
#' in each Delaunay triangle is forbidden for point generation).
#' @param DTmesh Delaunay triangulation of \code{Yp},
#' default is \code{NULL},
#' which is computed via \code{\link[interp]{tri.mesh}} function
#' in \code{interp} package.
#' \code{\link[interp]{tri.mesh}} function yields
#' the triangulation nodes with their neighbours, and
#' creates a triangulation object.
#' @param DTr Delaunay triangles based on \code{Yp},
#' default is \code{NULL},
#' which is computed via \code{\link[interp]{tri.mesh}} function
#' in \code{interp} package.
#' \code{\link[interp]{triangles}} function yields
#' a triangulation data structure from the triangulation object created
#' by \code{\link[interp]{tri.mesh}}.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title
#' for the plot of the point pattern}
#' \item{parameters}{Exclusion parameter, \code{delta},
#' of the Type I segregation pattern.
#' \code{delta} is in \eqn{(0,4/9)}
#' \eqn{\delta 100} \% area around each vertex in each Delaunay triangle
#' is forbidden for point generation.}
#' \item{ref.points}{The input set of points \code{Yp};
#' reference points, i.e.,
#' points from which generated points are segregated.}
#' \item{gen.points}{The output set of generated points segregated
#' from \code{Yp} points.}
#' \item{tri.Y}{Logical output, \code{TRUE},
#' if triangulation based on \code{Yp} points should be implemented.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points
#' and the number of reference (i.e., \code{Yp}) points.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points,
#' which are the \code{Yp} points}
#'
#' @seealso \code{\link{rseg.circular}}, \code{\link{rseg.std.tri}},
#' \code{\link{rsegII.std.tri}}, and \code{\link{rassoc.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-4;  #try also nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Yp<-cbind(runif(ny),runif(ny))
#' del<-.4
#'
#' Xdt<-rseg.multi.tri(nx,Yp,del)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' #or use
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
#' #Delaunay triangulation based on Y points
#' TRY<-interp::triangles(DTY)[,1:3];
#' Xp<-rseg.multi.tri(nx,Yp,del,DTY,TRY)$gen.points
#' #data under CSR in the convex hull of Ypoints
#'
#' Xlim<-range(Yp[,1])
#' Ylim<-range(Yp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' #plot of the data in the convex hull of Y points together with the Delaunay triangulation
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
#' #Delaunay triangulation based on Y points
#'
#' par(pty="s")
#' plot(Xp,main="Points from Type I Segregation \n in Multipe Triangles",
#' xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),
#' ylim=Ylim+yd*c(-.05,.05),type="n")
#' interp::plot.triSht(DTY, add=TRUE,
#' do.points=TRUE,col="blue")
#' points(Xp,pch=".",cex=3)
#' }
#'
#' @export rseg.multi.tri
rseg.multi.tri <- function(n,Yp,delta,DTmesh=NULL,DTr=NULL)
{
  Yp<-as.matrix(Yp)
  if (!is.numeric(Yp) || ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be numeric and of dimension kx2 with k>=3')}

  if (!is.point(delta,1) || delta<=0 || delta>=1)
  {stop('delta must be a scalar in (0,1)')}

  ny<-nrow(Yp)
  if (ny==3)
  {
    res<-rseg.tri(n,Yp,delta)
  } else
  {
    if (is.null(DTmesh))
    {DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    #Delaunay triangulation
    }

    if (is.null(DTr))
    {DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3);
    #the Delaunay triangles
    }
    DTr<-matrix(DTr,ncol=3)

    cnt<-1;ind<-0
    Xp<-NULL
    while (cnt <= n)
    {
      pnt<-as.vector(runif.multi.tri(1,Yp,DTmesh)$g)
     ind<-index.delaunay.tri(pnt,Yp,DTmesh)
      nodes<-as.numeric(DTr[ind,])
      tri<-Yp[nodes,]
      seg.tri<-seg.tri.support(delta,tri)
      if (in.triangle(pnt,seg.tri,boundary=TRUE)$in.tri)
      {
        Xp<-rbind(Xp,pnt);
        cnt<-cnt+1
      }
    }

        row.names(Xp)<-c()
    param<-delta
    names(param)<-"exclusion parameter"
    typ<-paste("Type I Segregation of ",n, " points from ",ny, " Y points with exclusion parameter delta = ",delta,sep="")

    npts<-c(n,ny)
    names(npts)<-c("nx","ny")

    txt<-"Type I Segregation of One Class from Another"
    main.txt<-paste("Type I Segregation in the Multi-Triangle Case \n with  Exclusion Parameter delta = ",delta,sep="")

    res<-list(
      type=typ,
      parameters=param,
      gen.points=Xp,
      #generated points segregated from Y points (vertices of Y triangles)
      ref.points=Yp,
      #reference points, i.e., points from which generated points are segregated (vertices of Y triangles)
      desc.pat=txt, #description of the pattern
      mtitle=main.txt,
      tri.Y=TRUE,
      num.points=npts,
      xlimit=range(Yp[,1]),
      ylimit=range(Yp[,2])
    )

    class(res)<-"Patterns"
    res$call <-match.call()
  }

  res
} #end of the function
#'

#################################################################

#' @title Generation of points associated (in a Type I fashion)
#' with a given set of points
#'
#' @description
#' An object of class \code{"Patterns"}.
#' Generates \code{n} points uniformly in the support
#' for Type I association in the convex hull of set of points, \code{Yp}.
#' \code{delta} is the parameter of association
#' (that is, only \eqn{\delta 100} \% area around each vertex in each Delaunay
#' triangle is allowed for point generation).
#'
#' \code{delta} corresponds to \code{eps}
#' in the standard equilateral triangle
#' \eqn{T_e} as \eqn{delta=4eps^2/3}
#' (see \code{rseg.std.tri}  function).
#'
#' If \code{Yp} consists only of 3 points,
#' then the function behaves like the
#' function \code{\link{rassoc.tri}}.
#'
#' \code{DTmesh} must be the Delaunay triangulation of \code{Yp}
#' and \code{DTr} must be the corresponding Delaunay triangles
#' (both \code{DTmesh} and \code{DTr} are \code{NULL} by default).
#' If \code{NULL}, \code{DTmesh} is computed via
#' \code{\link[interp]{tri.mesh}} and \code{DTr} is computed via
#' \code{\link[interp]{triangles}} function in \code{interp} package.
#'
#' \code{\link[interp]{tri.mesh}} function yields
#' the triangulation nodes with their neighbours,
#' and creates a triangulation object,
#' and \code{\link[interp]{triangles}} function yields
#' a triangulation data structure from the triangulation object created
#' by \code{\link[interp]{tri.mesh}}
#' (the first three columns are the vertex indices of the Delaunay triangles).
#'
#' See (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the association pattern.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param n A positive integer
#' representing the number of points to be generated.
#' @param Yp A set of 2D points
#' from which Delaunay triangulation is constructed.
#' @param delta A positive real number in \eqn{(0,4/9)}.
#' \code{delta} is the parameter of association (that is, only
#' \eqn{\delta 100} \% area around each vertex
#' in each Delaunay triangle is allowed for point generation).
#' @param DTmesh Delaunay triangulation of \code{Yp}, default is \code{NULL},
#' which is computed via \code{\link[interp]{tri.mesh}} function
#' in \code{interp} package. \code{\link[interp]{tri.mesh}} function yields
#' the triangulation nodes with their neighbours, and
#' creates a triangulation object.
#' @param DTr Delaunay triangles based on \code{Yp}, default is \code{NULL},
#' which is computed via \code{\link[interp]{tri.mesh}} function
#' in \code{interp} package. \code{\link[interp]{triangles}} function yields
#' a triangulation data structure from the triangulation object created
#' by \code{\link[interp]{tri.mesh}}.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{parameters}{Attraction parameter, \code{delta},
#' of the Type I association pattern.
#' \code{delta} is in \eqn{(0,4/9)}
#' only \eqn{\delta 100} \% of the area around each vertex in each Delaunay triangle
#' is allowed for point generation.}
#' \item{ref.points}{The input set of points \code{Yp};
#' reference points, i.e., points with which generated points are associated.}
#' \item{gen.points}{The output set of generated points
#' associated with \code{Yp} points.}
#' \item{tri.Y}{Logical output,
#' \code{TRUE} if triangulation based on
#' \code{Yp} points should be implemented.}
#' \item{desc.pat}{Description of the point pattern}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points
#' and the number of reference (i.e., \code{Yp}) points.}
#' \item{xlimit,ylimit}{The ranges of the \eqn{x}-
#' and \eqn{y}-coordinates of the reference points, which are the
#' \code{Yp} points}
#'
#' @seealso \code{\link{rassoc.circular}}, \code{\link{rassoc.std.tri}},
#' \code{\link{rassocII.std.tri}}, and \code{\link{rseg.multi.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Yp<-cbind(runif(ny),runif(ny))
#' del<-.4
#'
#' Xdt<-rassoc.multi.tri(nx,Yp,del)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' #or use
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
#' #Delaunay triangulation based on Y points
#' TRY<-interp::triangles(DTY)[,1:3];
#' Xp<-rassoc.multi.tri(nx,Yp,del,DTY,TRY)$g
#' #data under CSR in the convex hull of Ypoints
#'
#' Xlim<-range(Yp[,1])
#' Ylim<-range(Yp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' #plot of the data in the convex hull of Y points together with the Delaunay triangulation
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
#' #Delaunay triangulation based on Y points
#'
#' plot(Xp,main="Points from Type I Association \n in Multipe Triangles",
#' xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),
#' ylim=Ylim+yd*c(-.05,.05),type="n")
#' interp::plot.triSht(DTY, add=TRUE,
#' do.points=TRUE,col="blue")
#' points(Xp,pch=".",cex=3)
#' }
#'
#' @export rassoc.multi.tri
rassoc.multi.tri <- function(n,Yp,delta,DTmesh=NULL,DTr=NULL)
{
  Yp<-as.matrix(Yp)
  if (!is.numeric(Yp) || ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be numeric and of dimension kx2 with k>=3')}

  if (!is.point(delta,1) || delta<=0 || delta>=1)
  {stop('delta must be a scalar in (0,1)')}

  ny<-nrow(Yp)
  if (ny==3)
  {
    res<-rassoc.tri(n,Yp,delta)
  } else
  {
    if (is.null(DTmesh))
    {DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    #Delaunay triangulation
    }

    if (is.null(DTr))
    {DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3);
    #the Delaunay triangles
    }

    DTr<-matrix(DTr,ncol=3)

    cnt<-1;ind<-0
    Xp<-NULL
    while (cnt <= n)
    {
      pnt<-as.vector(runif.multi.tri(1,Yp,DTmesh)$g)
     ind<-index.delaunay.tri(pnt,Yp,DTmesh)
      nodes<-as.numeric(DTr[ind,])
      tri<-Yp[nodes,]
      seg.tri<-seg.tri.support(delta,tri)
      if (!in.triangle(pnt,seg.tri)$in.tri)
      {
        Xp<-rbind(Xp,pnt);
        cnt<-cnt+1
      }
    }

    row.names(Xp)<-c()
    param<-delta
    names(param)<-"attraction parameter"
    typ<-paste("Type I Association of ",n, " points with ",ny, " Y points with attraction parameter delta = ",delta,sep="")

    ny<-nrow(Yp)
    npts<-c(n,ny)
    names(npts)<-c("nx","ny")

    txt<-"Type I Association of One Class with Y points"
    main.txt<-paste("Type I Association in the Multi-Triangle Case \n with Attraction Parameter delta = ",delta,sep="")

    res<-list(
      type=typ,
      parameters=param,
      gen.points=Xp,
      #generated points associated Y points (vertices of Y triangles)
      ref.points=Yp,
      #reference points, i.e., points with which generated points are associated (vertices of Y triangles)
      desc.pat=txt, #description of the pattern
      mtitle=main.txt,
      tri.Y=TRUE,
      num.points=npts,
      xlimit=range(Yp[,1]),
      ylimit=range(Yp[,2])
    )

    class(res)<-"Patterns"
    res$call <-match.call()
  }

  res
} #end of the function
#'

#################################################################
#Functions in R^3
#################################################################

#' @title Generation of Uniform Points
#' in the Standard Regular Tetrahedron \eqn{T_h}
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly in the standard regular tetrahedron
#' \eqn{T_h=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))}.
#'
#' @param n A positive integer
#' representing the number of uniform points to be generated in the
#' standard regular tetrahedron \eqn{T_h}.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title for the plot of the point pattern}
#' \item{tess.points}{The vertices of the support region
#' of the uniformly generated points, it is the
#' standard regular tetrahedron \eqn{T_h} for this function}
#' \item{gen.points}{The output set of generated points uniformly
#' in the standard regular tetrahedron \eqn{T_h}.}
#' \item{out.region}{The outer region which contains the support region,
#' \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern
#' from which points are to be generated}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number
#' of vertices of the support points (here it is 4).}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}}
#' \item{xlimit,ylimit,zlimit}{The ranges of the \eqn{x}-, \eqn{y}-,
#' and \eqn{z}-coordinates of the support, \eqn{T_h}}
#'
#' @seealso \code{\link{runif.tetra}}, \code{\link{runif.tri}},
#' and \code{\link{runif.multi.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' tetra<-rbind(A,B,C,D)
#' n<-100
#'
#' set.seed(1)
#' Xdt<-runif.std.tetra(n)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-runif.std.tetra(n)$g
#'
#' Xlim<-range(tetra[,1])
#' Ylim<-range(tetra[,2])
#' Zlim<-range(tetra[,3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3],
#' phi =20,theta=15, bty = "g", pch = 20, cex = 1,
#' ticktype = "detailed",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),
#' zlim=Zlim+zd*c(-.05,.05))
#' #add the vertices of the tetrahedron
#' plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],
#' add=TRUE,lwd=2)
#'
#' plot3D::text3D(tetra[,1]+c(.05,0,0,0),tetra[,2],tetra[,3],
#' labels=c("A","B","C","D"), add=TRUE)
#' }
#'
#' \dontrun{
#' #need to install scatterplot3d package and call "library(scatterplot3d)"
#' s3d<-scatterplot3d(Xp, highlight.3d=TRUE,xlab="x",
#' ylab="y",zlab="z",
#'               col.axis="blue", col.grid="lightblue",
#'                main="3D Scatterplot of the data", pch=20)
#' s3d$points3d(tetra,pch=20,col="blue")
#' }
#'
#' @export runif.std.tetra
runif.std.tetra <- function(n)
{
  X <-matrix(0,n,3);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,0,1); y<-runif(1,0,.866025404); z<-runif(1,0,.816496581)
    if (y < 1.732050808*x & y < 1.732050808-1.732050808*x &
        z < -2.449489743*x+2.449489743-1.414213562*y &
        z < 2.449489743*x-1.414213562*y &
        z < 2.828427124*y)
    {X[i,]<-c(x,y,z);
    ct<-1;
    }
  }
  }
  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0);
  D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  Th<-rbind(A,B,C,D)  #std regular tetrahedron

  typ<-"Uniform Distribution in the Standard Regular Tetrahedron with Vertices A=(0,0,0), B=(1,0,0), C=(1/2,sqrt(3)/2,0), and D=(1/2,sqrt(3)/6,sqrt(6)/3)"
  main.txt<-"Uniform Points in the \n Standard Regular Tetrahedron"
  txt<-paste(n, " uniform points in the standard regular tetrahedron")

  ny<-4
  npts<-c(n,ny)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], "is the number of Uniform points \n",names(npts)[2], "is the number of vertices of the support region")

  res<-list(
    type=typ,
    mtitle=main.txt,
    tess.points=Th,
    #tessellation points whose convex hull constitutes the support of the uniform points
    gen.points=X,
    #uniformly generated points in the std regular tetrahedron
    out.region=NULL, #outer region for Te
    desc.pat=txt, #description of the pattern
    num.points=npts, txt4pnts=txt4pnts,
    xlimit=range(Th[,1]),
    ylimit=range(Th[,2]),
    zlimit=range(Th[,3])
  )

  class(res)<-"Uniform"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Generation of Uniform Points in a tetrahedron
#'
#' @description
#' An object of class \code{"Uniform"}.
#' Generates \code{n} points uniformly
#' in the general tetrahedron \code{th}
#' whose vertices are stacked row-wise.
#'
#' @param n A positive integer
#' representing the number of uniform points
#' to be generated in the tetrahedron.
#' @param th A \eqn{4 \times 3} matrix with each row
#' representing a vertex of the tetrahedron.
#'
#' @return A \code{list} with the elements
#' \item{type}{The type of the pattern
#' from which points are to be generated}
#' \item{mtitle}{The \code{"main"} title
#' for the plot of the point pattern}
#' \item{tess.points}{The vertices of the support of
#' the uniformly generated points,
#' it is the tetrahedron' \code{th} for this function}
#' \item{gen.points}{The output set of generated points uniformly
#' in the tetrahedron, \code{th}.}
#' \item{out.region}{The outer region
#' which contains the support region,
#' \code{NULL} for this function.}
#' \item{desc.pat}{Description of the point pattern
#' from which points are to be generated}
#' \item{num.points}{The \code{vector} of two numbers,
#' which are the number of generated points and the number
#' of vertices of the support points (here it is 4).}
#' \item{txt4pnts}{Description of the two numbers in \code{num.points}}
#' \item{xlimit, ylimit, zlimit}{The ranges of the \eqn{x}-, \eqn{y}-,
#' and \eqn{z}-coordinates of the support, \code{th}}
#'
#' @seealso \code{\link{runif.std.tetra}} and \code{\link{runif.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-sample(1:12,3); B<-sample(1:12,3);
#' C<-sample(1:12,3); D<-sample(1:12,3)
#' tetra<-rbind(A,B,C,D)
#'
#' n<-100
#'
#' set.seed(1)
#' Xdt<-runif.tetra(n,tetra)
#' Xdt
#' summary(Xdt)
#' plot(Xdt)
#'
#' Xp<-Xdt$g
#'
#' Xlim<-range(tetra[,1],Xp[,1])
#' Ylim<-range(tetra[,2],Xp[,2])
#' Zlim<-range(tetra[,3],Xp[,3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3],
#' theta =225, phi = 30, bty = "g",
#' main="Uniform Points in a Tetrahedron",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),
#' zlim=Zlim+zd*c(-.05,.05),
#'           pch = 20, cex = 1, ticktype = "detailed")
#' #add the vertices of the tetrahedron
#' plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],
#' add=TRUE,lwd=2)
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3],
#' labels=c("A","B","C","D"), add=TRUE)
#'
#' #need to install scatterplot3d package and call "library(scatterplot3d)"
#' s3d<-scatterplot3d(Xp, highlight.3d=TRUE,
#' xlab="x",ylab="y",zlab="z",
#'                     col.axis="blue", col.grid="lightblue",
#'                     main="3D Scatterplot of the data", pch=20)
#'  s3d$points3d(tetra,pch=20,col="blue")
#'  }
#'
#' @export runif.tetra
runif.tetra <- function(n,th)
{
  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  X <-matrix(0,n,3);
  for (i in 1:n)
  {ct<-0;
  while (ct==0)
  {
    x<-runif(1,min(th[,1]),max(th[,1]));
    y<-runif(1,min(th[,2]),max(th[,2]));
    z<-runif(1,min(th[,3]),max(th[,3]))
    p<-c(x,y,z)
    if (in.tetrahedron(p,th)$in.tetra)
    {X[i,]<-p;
    ct<-1;
    }
  }
  }

  Avec= paste(round(th[1,],2), collapse=",");
  Bvec= paste(round(th[2,],2), collapse=",");
  Cvec= paste(round(th[3,],2), collapse=",");
  Dvec= paste(round(th[4,],2), collapse=",");
  main.txt<-paste("Uniform Points in the Tetrahedron with Vertices \n (",Avec,"), (",Bvec,"), (",Cvec,") and (",Dvec,")",sep="")
  txt<-paste(n, " uniform points in the tetrahedron with vertices (",Avec,"), (",Bvec,"), (",Cvec,") and (",Dvec,")",sep="")
  typ<-paste("Uniform Distribution in the Tetrahedron with Vertices (",Avec,"), (",Bvec,"), (",Cvec,") and (",Dvec,")",sep="")

  ny<-4
  npts<-c(n,ny)
  names(npts)<-c("nx","ny")
  txt4pnts<-paste(names(npts)[1], "is the number of Uniform points \n",names(npts)[2], "is the number of vertices of the support region")

  res<-list(
    type=typ,
    mtitle=main.txt,
    gen.points=X,
    #uniformly generated points in the std regular tetrahedron
    tess.points=th,
    #tessellation points whose convex hull constitutes the support of the uniform points
    desc.pat=txt, #description of the pattern
    out.region=NULL, #outer region for th
    num.points=npts,  txt4pnts=txt4pnts,
    xlimit=range(th[,1]),
    ylimit=range(th[,2]),
    zlimit=range(th[,3])
  )

  class(res)<-"Uniform"
  res$call <-match.call()
  res
} #end of the function
#'
