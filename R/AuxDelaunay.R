#AuxDelaunay.R
#Contains the auxiliary functions used in PCD calculations, such as edge and vertex regions
#####AUXILIARY FUNCTIONS#################

#' @import graphics
#' @import stats
#' @import plot3D
#' @import interp
#' @importFrom Rdpack reprompt

#################################################################
#Auxiliary Function for Triangles
#################################################################

# funsAB2CMTe
#'
#' @title The lines joining two vertices to the center of mass in standard equilateral triangle
#'
#' @description Two functions, \code{lA_CM.Te} and \code{lB_CM.Te} of class \code{"TriLines"}.
#' Returns the \code{equation, slope, intercept}, and \eqn{y}-coordinates of the lines joining \eqn{A} and \eqn{CM} and
#' also \eqn{B} and \eqn{CM}.
#'
#' \code{lA_CM.Te} is the line joining \eqn{A} to the center of mass, \eqn{CM}, and
#'
#' \code{lB_CM.Te} is the line joining \eqn{B} to the center of mass, \eqn{CM}, in the standard equilateral triangle \eqn{T_e=(A,B,C)}
#' with \eqn{A=(0,0)}, \eqn{B=(1,0)}, \eqn{C=(1/2,\sqrt{3}/2)}; \eqn{x}-coordinates are provided in \code{vector} \code{x}.
#'
#' @param x A single scalar or a \code{vector} of scalars which is the argument of the functions
#' \code{lA_CM.Te} and \code{lB_CM.Te}.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Longer description of the line.}
#' \item{txt2}{Shorter description of the line (to be inserted over the line in the plot).}
#' \item{mtitle}{The \code{"main"} title for the plot of the line.}
#' \item{cent}{The center chosen inside the standard equilateral triangle.}
#' \item{cent.name}{The name of the center inside the standard equilateral triangle.
#' It is \code{"CM"} for these two functions.}
#' \item{tri}{The triangle (it is the standard equilateral triangle for this function).}
#' \item{x}{The input vector, can be a scalar or a \code{vector} of scalars,
#' which constitute the \eqn{x}-coordinates of the point(s) of interest on the line.}
#' \item{y}{The output vector, will be a scalar if \code{x} is a scalar or a \code{vector} of scalars
#' if \code{x} is a \code{vector} of scalar, constitutes the \eqn{y}-coordinates of the point(s) of interest on the line.}
#' \item{slope}{Slope of the line.}
#' \item{intercept}{Intercept of the line.}
#' \item{equation}{Equation of the line.}
#'
#' @name funsAB2CMTe
NULL
#'
#' @seealso \code{\link{lA_M.Te}}, \code{\link{lB_M.Te}}, and \code{\link{lC_M.Te}}
#'
#' @rdname funsAB2CMTe
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for lA_CM.Te
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.1)  #try also by=.01
#'
#' lnACM<-lA_CM.Te(x)
#' lnACM
#' summary(lnACM)
#' plot(lnACM)
#'
#' CM<-(A+B+C)/3;
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' L<-Te; R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Te,CM,D1,D2,D3,c(.25,lA_CM.Te(.25)$y),c(.75,lB_CM.Te(.75)$y))
#' xc<-txt[,1]+c(-.02,.02,.02,.05,.05,-.03,.0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,.02,0,.02,-.04,0,0)
#' txt.str<-c("A","B","C","CM","D1","D2","D3","lA_CM.Te(x)","lB_CM.Te(x)")
#' text(xc,yc,txt.str)
#'
#' lA_CM.Te(.25)$y
#' }
#'
#' @export
lA_CM.Te <- function(x)
{
  dname <-deparse(substitute(x))

  if (!is.vector(x))
  {stop('the argument must be a vector')}

  sl<-.5773502694
  int<-0
  ln<-int+sl*x

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2); #the vertices of the standard equilateral triangle
  Tr<-rbind(A,B,C)
  M<-(A+B+C)/3;

  names(sl)<-"slope"
  names(int)<-"intercept"

  cname <-"CM"

  txt1<-"Line Passing through A and the Center of Mass, CM, in the Standard Equilateral Triangle Te=ABC\n with A=(0,0);  B=(1,0); C=(1/2,sqrt(3)/2)"
  main.text<-"Line Passing through A and the Center of Mass CM\n in the Standard Equilateral Triangle"
  txt2<-paste("lA_CM.Te(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2, mtitle=main.text,
    cent=M, cent.name=cname,
    tri=Tr,
    x=x,
    y=ln,
    slope=sl,
    intercept=int,
    equation=ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                                                             ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'
#' @rdname funsAB2CMTe
#'
#' @examples
#' \dontrun{
#' #Examples for lB_CM.Te
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.1)  #try also by=.01
#'
#' lnBCM<-lB_CM.Te(x)
#' lnBCM
#' summary(lnBCM)
#' plot(lnBCM,xlab="x",ylab="y")
#'
#' lB_CM.Te(.25)$y
#' }
#'
#' @export
lB_CM.Te <- function(x)
{
  dname <-deparse(substitute(x))

  if (!is.vector(x))
  {stop('the argument must be a vector')}

  sl<--.5773502694
  int<-.5773502694
  ln<-int+sl*x

  names(sl)<-"slope"
  names(int)<-"intercept"

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2); #the vertices of the standard equilateral triangle
  Tr<-rbind(A,B,C)
  M<-(A+B+C)/3;

  cname <-"CM"

  txt1<-"Line Passing through B and the Center of Mass, CM, in the Standard Equilateral Triangle Te=ABC\n with A=(0,0);  B=(1,0); C=(1/2,sqrt(3)/2)"
  main.text<-"Line Passing through B and the Center of Mass CM\n in the Standard Equilateral Triangle"
  txt2<-paste("lB_CM.Te(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2, mtitle=main.text,
    cent=M, cent.name=cname,
    tri=Tr,
    x=x,
    y=ln,
    slope=sl,
    intercept=int,
    equation=ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                                                             ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'
#################################################################

# funsAB2MTe
#'
#' @title The lines joining the three vertices of the standard equilateral triangle to a center, \code{M}, of it
#'
#' @description
#' Three functions, \code{lA_M.Te}, \code{lB_M.Te} and \code{lC_M.Te} of class \code{"TriLines"}.
#' Returns the \code{equation, slope, intercept}, and \eqn{y}-coordinates of the lines joining \eqn{A} and \code{M},
#' \eqn{B} and \code{M}, and also \eqn{C} and \code{M}.
#'
#' \code{lA_M.Te} is the line joining \eqn{A} to the center, \code{M}, \code{lB_M.Te} is the line joining \eqn{B} to \code{M},
#' and \code{lC_M.Te} is the line joining C to \code{M}, in the standard equilateral triangle \eqn{T_e=(A,B,C)}
#' with \eqn{A=(0,0)}, \eqn{B=(1,0)}, \eqn{C=(1/2,\sqrt{3}/2)}; \eqn{x}-coordinates are provided in \code{vector} \code{x}
#'
#' @param x A single scalar or a \code{vector} of scalars.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Longer description of the line.}
#' \item{txt2}{Shorter description of the line (to be inserted over the line in the plot).}
#' \item{mtitle}{The \code{"main"} title for the plot of the line.}
#' \item{cent}{The center chosen inside the standard equilateral triangle.}
#' \item{cent.name}{The name of the center inside the standard equilateral triangle.}
#' \item{tri}{The triangle (it is the standard equilateral triangle for this function).}
#' \item{x}{The input vector, can be a scalar or a \code{vector} of scalars,
#' which constitute the \eqn{x}-coordinates of the point(s) of interest on the line.}
#' \item{y}{The output vector, will be a scalar if \code{x} is a scalar or a \code{vector} of scalars if \code{x} is
#' a \code{vector} of scalar, constitutes the \eqn{y}-coordinates of the point(s) of interest on the line.}
#' \item{slope}{Slope of the line.}
#' \item{intercept}{Intercept of the line.}
#' \item{equation}{Equation of the line.}
#'
#' @name funsAB2MTe
NULL
#'
#' @seealso \code{\link{lA_CM.Te}} and \code{\link{lB_CM.Te}}
#'
#' @rdname funsAB2MTe
#'
#' @examples
#' \dontrun{
#' #Examples for lA_M.Te
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#'
#' M<-c(.65,.2)  #try also M<-c(1,1,1)
#'
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.1)  #try also by=.01
#'
#' lnAM<-lA_M.Te(x,M)
#' lnAM
#' summary(lnAM)
#' plot(lnAM)
#'
#' Ds<-cent2edges.tri(Te,M)
#' #finds the projections from a point M=(m1,m2) to the edges on the
#' #extension of the lines joining M to the vertices in the triangle Te
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,pch=".",xlab="",ylab="",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' L<-Te; R<-rbind(M,M,M)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' L<-Ds; R<-rbind(M,M,M)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=3,col=2)
#'
#' txt<-rbind(Te,M,Ds,c(.25,lA_M.Te(.25,M)$y),c(.4,lB_M.Te(.4,M)$y),
#' c(.60,lC_M.Te(.60,M)$y))
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.04,-.03,.0,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,.05,.02,.03,-.03,0,0,0)
#' txt.str<-c("A","B","C","M","D1","D2","D3","lA_M.Te(x)","lB_M.Te(x)","lC_M.Te(x)")
#' text(xc,yc,txt.str)
#'
#' lA_M.Te(.25,M)
#' }
#'
#' @export
lA_M.Te <- function(x,M)
{
  dname <-deparse(substitute(x))
  cname <-deparse(substitute(M))

  if (!is.point(x,length(x)))
  {stop('x must be a numeric vector')}

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {
    M<-bary2cart(M,Te)
    cname <-"M"
  }
  if (!in.triangle(M,Te,boundary = FALSE)$in.tri)
  {stop('M must be in the interior of the standard equilateral triangle')}

  m1<-M[1]
  m2<-M[2]

  sl<-m2/m1
 int<-m2+(0-m1)*m2/m1
  ln<-int+sl*x

  names(sl)<-"slope"
  names(int)<-"intercept"

  txt1<-"Line Passing through A and Center, M, in the Standard Equilateral Triangle Te=ABC\n with A=(0,0);  B=(1,0); C=(1/2,sqrt(3)/2)"
  main.text<-"Line Passing through A and the Center M\n in the Standard Equilateral Triangle"
  txt2<-paste("lA_M.Te(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2, mtitle=main.text,
    cent=M, cent.name=cname,
    tri=Te,
    x=x,
    y=ln,
    slope=sl,
   intercept=int,
    equation=ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                                                             ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'
#' @rdname funsAB2MTe
#'
#' @examples
#' \dontrun{
#' #Examples for lB_M.Te
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#'
#' M<-c(.65,.2)  #try also M<-c(1,1,1)
#'
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.5)  #try also by=.1
#'
#' lnBM<-lB_M.Te(x,M)
#' lnBM
#' summary(lnBM)
#' plot(lnBM)
#' }
#'
#' @export
lB_M.Te <- function(x,M)
{
  dname <-deparse(substitute(x))
  cname <-deparse(substitute(M))

  if (!is.point(x,length(x)))
  {stop('x must be a numeric vector')}

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)
  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {
    M<-bary2cart(M,Te)
    cname <-"M"
  }
  if (!in.triangle(M,Te,boundary = FALSE)$in.tri)
  {stop('M must be in the interior of the standard equilateral triangle')}

  m1<-M[1]
  m2<-M[2]

  sl<-m2/(m1-1)
 int<-m2+(0-m1)*m2/(m1-1)
  ln<-int+sl*x

  names(sl)<-"slope"
  names(int)<-"intercept"

  txt1<-"Line Passing through B and the Center, M, in the Standard Equilateral Triangle Te=ABC\n with A=(0,0);  B=(1,0); C=(1/2,sqrt(3)/2)"
  main.text<-"Line Passing through B and the Center M\n in the Standard Equilateral Triangle"
  txt2<-paste("lB_M.Te(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2, mtitle=main.text,
    cent=M, cent.name=cname,
    tri=Te,
    x=x,
    y=ln,
    slope=sl,
   intercept=int,
    equation=ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                                                             ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'
#' @rdname funsAB2MTe
#'
#' @examples
#' \dontrun{
#' #Examples for lC_M.Te
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#'
#' M<-c(.65,.2)  #try also M<-c(1,1,1)
#'
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.5)  #try also by=.1
#'
#' lnCM<-lC_M.Te(x,M)
#' lnCM
#' summary(lnCM)
#' plot(lnCM)
#' }
#'
#' @export
lC_M.Te <- function(x,M)
{
  dname <-deparse(substitute(x))
  cname <-deparse(substitute(M))

  if (!is.vector(x))
  {stop('x must be a vector')}

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)
  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}


  if (dimension(M)==3)
  {
    M<-bary2cart(M,Te)
    cname <-"M"
  }
  if (!in.triangle(M,Te,boundary = FALSE)$in.tri)
  {stop('M must be in the interior of the standard equilateral triangle')}

  m1<-M[1];
  m2<-M[2];

  sl<-(m2*2-sqrt(3))/(2*m1-1)
 int<-(m2*(2*0-1)-sqrt(3)*(0-m1))/(2*m1-1)
  ln<-int+sl*x

  names(sl)<-"slope"
  names(int)<-"intercept"

  txt1<-"Line Passing through C and the Center, M, in the Standard Equilateral Triangle Te=ABC\n with A=(0,0);  B=(1,0); C=(1/2,sqrt(3)/2)"
  main.text<-"Line Passing through C and the Center M\n in the Standard Equilateral Triangle"
  txt2<-paste("lC_M.Te(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2, mtitle=main.text,
    cent=M, cent.name=cname,
    tri=Te,
    x=x,
    y=ln,
    slope=sl,
   intercept=int,
    equation=ifelse(sl==0 & int==0,"y=0",ifelse(sl!=0 & int==0,paste("y=",sl,"x",sep=""),ifelse(sl==0 & int!=0,paste("y=",int,sep=""),ifelse(sl==1,ifelse(sign(int)<0,paste("y=x",int,sep=""),paste("y=x+",int,sep="")),
                                                                                                                                             ifelse(sign(int)<0,paste("y=",sl,"x",int,sep=""),paste("y=",sl,"x+",int,sep="")))) ))
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title Circumcenter of a standard basic triangle form
#'
#' @description Returns the circumcenter of a standard basic triangle form \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} given \eqn{c_1}, \eqn{c_2}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle form by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See (\insertCite{weisstein-tri-centers,ceyhan:comp-geo-2010;textual}{pcds}) for triangle centers and
#' (\insertCite{ceyhan:arc-density-PE,ceyhan:arc-density-CS,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}) for the
#' standard basic triangle form.
#'
#' @param c1,c2 Positive real numbers representing the top vertex in standard basic triangle form \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' @return circumcenter of the standard basic triangle form \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} given \eqn{c_1}, \eqn{c_2} as the arguments of the
#' function
#'
#' @seealso \code{\link{circ.cent.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);  #the vertices of the standard basic triangle form Tb
#' Tb<-rbind(A,B,C)
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' CC
#'
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tb[,1])
#' Ylim<-range(Tb[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' par(pty="s")
#' plot(A,pch=".",asp=1,xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' points(rbind(CC))
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tb,CC,D1,D2,D3)
#' xc<-txt[,1]+c(-.03,.04,.03,.06,.06,-.03,0)
#' yc<-txt[,2]+c(.02,.02,.03,-.03,.02,.04,-.03)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' #for an obtuse triangle
#' c1<-.4; c2<-.3;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);  #the vertices of the standard basic triangle form Tb
#' Tb<-rbind(A,B,C)
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' CC
#'
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tb[,1],CC[1])
#' Ylim<-range(Tb[,2],CC[2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' par(pty="s")
#' plot(A,pch=".",asp=1,xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' points(rbind(CC))
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tb,CC,D1,D2,D3)
#' xc<-txt[,1]+c(-.03,.03,.03,.07,.07,-.05,0)
#' yc<-txt[,2]+c(.02,.02,.04,-.03,.03,.04,.06)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export circ.cent.bas.tri
circ.cent.bas.tri <- function(c1,c2)
{
  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('both arguments must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  c(1/2, 1/2*(c2^2-c1+c1^2)/c2)
} #end of the function
#'

#################################################################

#' @title Circumcenter of a general triangle
#'
#' @description Returns the circumcenter a given triangle, \code{tri}, with vertices stacked row-wise.
#' See (\insertCite{weisstein-tri-centers,ceyhan:comp-geo-2010;textual}{pcds}) for triangle centers.
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return circumcenter of the triangle \code{tri}
#'
#' @seealso \code{\link{circ.cent.bas.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);  #the vertices of the triangle Tr
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#' CC
#'
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],CC[1])
#' Ylim<-range(Tr[,2],CC[2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,asp=1,pch=".",xlab="",ylab="",main="Circumcenter of a triangle",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(rbind(CC))
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tr,CC,Ds)
#' xc<-txt[,1]+c(-.08,.08,.08,.12,-.09,-.1,-.09)
#' yc<-txt[,2]+c(.02,-.02,.03,-.06,.02,.06,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);  #the vertices of the equilateral triangle Te
#' circ.cent.tri(Te)  #the circumcenter
#'
#' A<-c(0,0); B<-c(0,1); C<-c(2,0);
#' Tr<-rbind(A,B,C);  #the vertices of the triangle T
#' circ.cent.tri(Tr)  #the circumcenter
#' }
#'
#' @export circ.cent.tri
circ.cent.tri <- function(tri)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('the argument must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,];
  a1<-A[1]; a2<-A[2];
  b1<-B[1]; b2<-B[2];
  c1<-C[1]; c2<-C[2];

  Dx<--det(matrix(cbind(c(sum(A^2),sum(B^2),sum(C^2)),tri[,2],vec1),ncol=3))
  Dy<-det(matrix(cbind(c(sum(A^2),sum(B^2),sum(C^2)),tri[,1],vec1),ncol=3))

  cc1<--Dx/(2*D0); cc2<--Dy/(2*D0);
  c(cc1,cc2)
} #end of the function
#'

#################################################################

#' @title Projections of a point inside a triangle to its edges
#'
#' @description Returns the projections from a general center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
#' barycentric coordinates in the interior of a triangle to the edges on the extension of the lines joining
#' \code{M} to the vertices (see the examples for an illustration).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}.
#'
#' @return Three projection points (stacked row-wise) from a general center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of a triangle to the edges on
#' the extension of the lines joining \code{M} to the vertices;
#' row \eqn{i} is the projection  point into edge \eqn{i}, for \eqn{i=1,2,3}.
#'
#' @seealso \code{\link{cent2edges.bas.tri}} and \code{\link{cent.ndPE2edges}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' Ds<-cent2edges.tri(Tr,M)  #try also cent2edges.tri(Tr,M=c(1,1))
#' Ds
#'
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Tr,pch=".",xlab="",ylab="",main="Projection of Center M on the edges of a triangle",axes=TRUE,
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]
#' yc<-Tr[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(-.02,.04,-.04,-.02)
#' yc<-txt[,2]+c(-.02,.04,.04,-.06)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cent2edges.tri
cent2edges.tri <- function(tri,M)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];
  a1<-y1[1]; a2<-y1[2]; b1<-y2[1]; b2<-y2[2]; c1<-y3[1]; c2<-y3[2];

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  m1<-M[1]; m2<-M[2]

  M11<--(a1*b1*c2-a1*b1*m2-a1*b2*c1+a1*c1*m2+a2*b1*m1-a2*c1*m1-b1*c2*m1+b2*c1*m1)/(b2*a1-c2*a1-a2*b1+a2*c1+m2*b1-b2*m1-m2*c1+c2*m1);
  M12<-(a1*b2*m2-a1*c2*m2-a2*b1*c2+a2*b2*c1-a2*b2*m1+a2*c2*m1+b1*c2*m2-b2*c1*m2)/(b2*a1-c2*a1-a2*b1+a2*c1+m2*b1-b2*m1-m2*c1+c2*m1);
  D1<-c(M11,M12)
  M21 <-(a1*b1*c2-a1*b1*m2+a1*b2*m1-a1*c2*m1-a2*b1*c1+a2*c1*m1+b1*c1*m2-b2*c1*m1)/(b2*a1-m2*a1-a2*b1+m1*a2+c2*b1-c1*b2+m2*c1-c2*m1);
  M22 <-(a1*b2*c2-a1*c2*m2-a2*b1*m2-a2*b2*c1+a2*b2*m1+a2*c1*m2+b1*c2*m2-b2*c2*m1)/(b2*a1-m2*a1-a2*b1+m1*a2+c2*b1-c1*b2+m2*c1-c2*m1);
  D2<-c(M21,M22)
  M31 <- (a1*b2*c1-a1*b2*m1-a1*c1*m2+a1*c2*m1-a2*b1*c1+a2*b1*m1+b1*c1*m2-b1*c2*m1)/(c2*a1-m2*a1-a2*c1+m1*a2-c2*b1+m2*b1+c1*b2-b2*m1);
  M32 <- (a1*b2*c2-a1*b2*m2-a2*b1*c2+a2*b1*m2-a2*c1*m2+a2*c2*m1+b2*c1*m2-b2*c2*m1)/(c2*a1-m2*a1-a2*c1+m1*a2-c2*b1+m2*b1+c1*b2-b2*m1);
  D3<-c(M31,M32)

  Ds<-rbind(D1,D2,D3)
  row.names(Ds)<-c()

  Ds
} #end of the function
#'

#################################################################

#' @title Projections of a point inside the standard basic triangle form to its edges
#'
#' @description Returns the projections from a general center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
#' barycentric coordinates in the interior of the standard basic triangle form \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' to the edges on the extension of the lines joining \code{M} to the vertices (see the examples for an illustration).
#' In the standard basic triangle form \eqn{T_b}, \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle form by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle form is useful for simulation studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle form
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle form.
#'
#' @return Three projection points (stacked row-wise) from a general center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of a standard basic triangle form to the edges on
#' the extension of the lines joining \code{M} to the vertices;
#' row \eqn{i} is the projection  point into edge \eqn{i}, for \eqn{i=1,2,3}.
#'
#' @seealso \code{\link{cent2edges.tri}} and \code{\link{cent.ndPE2edges}}
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
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' Ds<-cent2edges.bas.tri(c1,c2,M)
#' Ds
#'
#' Xlim<-range(Tb[,1])
#' Ylim<-range(Tb[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tb)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Tb,pch=".",xlab="",ylab="",axes=TRUE,
#' xlim=Xlim+xd*c(-.1,.1),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' L<-rbind(M,M,M); R<-Tb
#' segments(L[,1], L[,2], R[,1], R[,2], lty=3,col=2)
#'
#' xc<-Tb[,1]+c(-.04,.05,.04)
#' yc<-Tb[,2]+c(.02,.02,.03)
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(-.02,.03,-.03,0)
#' yc<-txt[,2]+c(-.02,.02,.02,-.03)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cent2edges.bas.tri
cent2edges.bas.tri <- function(c1,c2,M)
{
  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (in.triangle(M,Tb,boundary=FALSE)$in.tri==FALSE)
  {stop('center, M, is not in the interior of the standard basic triangle')}

  m1<-M[1]; m2<-M[2]

  M11<--c2/(m2*c1-c2*m1-m2)*m1;
  M12<--m2*c2/(m2*c1-c2*m1-m2);
  M1<-c(M11,M12)
  M21 <-m2/(m2*c1-c2*m1+c2)*c1;
  M22 <-m2*c2/(m2*c1-c2*m1+c2);
  M2<-c(M21,M22)
  M31 <- -(m2*c1-c2*m1)/(-m2+c2);
  M32 <- 0;
  M3<-c(M31,M32)
  prj.pnts<-rbind(M1,M2,M3)
  row.names(prj.pnts)<-c()
  prj.pnts
} #end of the function
#'

#################################################################

#' @title Projections of Centers for non-degenerate asymptotic distribution of domination number of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) to its edges
#'
#' @description Returns the projections from  center \code{cent} to the edges on the extension of the lines joining \code{cent} to the vertices
#' in the triangle, \code{tri}. Here M is one of the three centers which gives nondegenerate asymptotic distribution
#' of the domination number of PE-PCD for uniform data in \code{tri} for a given expansion parameter \code{r} in \eqn{(1,1.5]}.
#' The center label \code{cent} values \code{1,2,3} correspond to the vertices \eqn{M_1}, \eqn{M_2}, and \eqn{M_3} (i.e., row numbers in the
#' output of \code{cent.nondegPE(tri,r)}); default for \code{cent} is 1. \code{cent} becomes center of mass \eqn{CM} for \eqn{r=1.5}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]} for this function.
#' @param cent Index of the center (as \eqn{1,2,3} corresponding to \eqn{M_1,\,M_2,\,M_3}) which gives nondegenerate asymptotic
#' distribution of the domination number of PE-PCD for uniform data in \code{tri} for expansion parameter \code{r} in \eqn{(1,1.5]};
#' default \code{cent=1}.
#'
#' @return Three projection points (stacked row-wise) from one of the centers (as \eqn{1,2,3} corresponding to
#' \eqn{M_1,\,M_2,\,M_3}) which gives nondegenerate asymptotic distribution of the domination number of PE-PCD for
#' uniform data in \code{tri} for expansion parameter \code{r} in \eqn{(1,1.5]}.
#'
#' @seealso \code{\link{cent2edges.bas.tri}} and \code{\link{cent2edges.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' r<-1.35
#'
#' cent.ndPE2edges(Tr,r,cent=2)
#'
#' Ms<-cent.nondegPE(Tr,r)
#' M1=Ms[1,]
#'
#' Ds<-cent.ndPE2edges(Tr,r,cent=1)
#'
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",
#' main="Projections from a non-degeneracy center\n to the edges of the triangle",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Ms,pch=".",col=1)
#' polygon(Ms,lty=2)
#'
#' xc<-Tr[,1]+c(-.02,.03,.02)
#' yc<-Tr[,2]+c(-.02,.04,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' txt<-Ms
#' xc<-txt[,1]+c(-.02,.04,-.04)
#' yc<-txt[,2]+c(-.02,.04,.04)
#' txt.str<-c("M1","M2","M3")
#' text(xc,yc,txt.str)
#'
#' points(Ds,pch=4,col=2)
#' L<-rbind(M1,M1,M1); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2,lwd=2,col=4)
#' txt<-Ds
#' xc<-txt[,1]+c(-.02,.04,-.04)
#' yc<-txt[,2]+c(-.02,.04,.04)
#' txt.str<-c("D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' cent.ndPE2edges(Tr,r,cent=3)
#' #gives an error message if center index, cent, is different from 1, 2 or 3
#' cent.ndPE2edges(Tr,r=1.49,cent=2)
#' #gives an error message if r>1.5
#' }
#'
#' @export cent.ndPE2edges
cent.ndPE2edges <- function(tri,r,cent=1)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(r,1) || r<=1 || r>1.5)
  {stop('r must be a scalar in (1,1.5]')}

  if (cent!=1 & cent!=2 & cent!=3)
  {stop('center index, cent, must be 1, 2 or 3')}

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,];

  if (cent==3)
  {
    D2<-A+(2-r)*(C-A)
    D1<-B+(2-r)*(C-B)
    D3<-(A+B)/2
  }
  else
  {
    if (cent==1)
    {
      D3<-B+(2-r)*(A-B)
      D2<-C+(2-r)*(A-C)
      D1<-(B+C)/2
    }
    else
    {
      D3<-A+(2-r)*(B-A)
      D1<-C+(2-r)*(B-C)
      D2<-(A+C)/2
    }
  }
  Ds<-rbind(D1,D2,D3)
  row.names(Ds)<-c()

  Ds
} #end of the function
#'


#################################################################

#' @title Check whether a point is inside a triangle
#'
#' @description Checks if the point \code{p} lies in the triangle, \code{tri}, using the barycentric
#' coordinates, generally denoted as \eqn{(\alpha,\beta,\gamma)}.
#'
#' If all (normalized or non-normalized)
#' barycentric coordinates are positive then the point \code{p} is inside the triangle,
#' if all are nonnegative with one or more are zero, then \code{p} falls in the boundary. If some of the
#' barycentric coordinates are negative, then \code{p} falls outside the triangle.
#'
#' \code{boundary} is a logical argument (default=\code{TRUE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if the point, \code{p}, lies in the closure of the triangle (i.e., interior and boundary
#' combined); else, it checks if \code{p} lies in the interior of the triangle.
#'
#' @param p A 2D point to be checked whether it is inside the triangle or not.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param boundary A logical parameter (default=\code{TRUE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if the point, \code{p}, lies in the closure of the triangle (i.e., interior and boundary
#' combined); else, it checks if \code{p} lies in the interior of the triangle.
#'
#' @return A \code{list} with two elements
#' \item{in.tri}{A logical output, it is \code{TRUE}, if the point, \code{p}, is inside the triangle, \code{tri},
#' else it is \code{FALSE}.}
#' \item{barycentric}{The barycentric coordinates \eqn{(\alpha,\beta,\gamma)} of the point \code{p} with respect to
#' the triangle, \code{tri}.}
#'
#' @seealso \code{\link{in.tri.all}} and \code{\link[interp]{on.convex.hull}} from the \code{interp} package for documentation for \code{in.convex.hull}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2); p<-c(1.4,1.2)
#' Tr<-rbind(A,B,C)
#' in.triangle(p,Tr)
#'
#' p<-c(.4,-.2)
#' in.triangle(p,Tr)
#'
#' #for the vertex A
#' in.triangle(A,Tr)
#' in.triangle(A,Tr,boundary=FALSE)
#'
#' #for a point on the edge AB
#' D3<-(A+B)/2
#' in.triangle(D3,Tr)
#' in.triangle(D3,Tr,boundary=FALSE)
#'
#' #for a NA entry point
#' p<-c(NA,.2)
#' in.triangle(p,Tr)
#' }
#'
#' @export
in.triangle <- function(p,tri,boundary=TRUE)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  p1<-tri[1,]; p2<-tri[2,]; p3<-tri[3,];
  if (is.na(p[1]) || is.na(p[2]))
  {ind.tri<-FALSE; alpha<-beta<-gamma<-NA;
  } else
  {
    p.x<-p[1]; p.y<-p[2]
    p1.x<-p1[1]; p1.y<-p1[2]
    p2.x<-p2[1]; p2.y<-p2[2]
    p3.x<-p3[1]; p3.y<-p3[2]
    # barycentric coordinates
    alpha = ((p2.y - p3.y)*(p.x - p3.x) + (p3.x - p2.x)*(p.y - p3.y)) /
      ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y));
    beta = ((p3.y - p1.y)*(p.x - p3.x) + (p1.x - p3.x)*(p.y - p3.y)) /
      ((p2.y - p3.y)*(p1.x - p3.x) + (p3.x - p2.x)*(p1.y - p3.y));
    gamma = 1 - alpha - beta;
    #if all of alpha, beta, and gamma are greater than 0,
    #then the point p lies within the interior of the triangle made of points p1, p2, and p3.

    if (boundary==TRUE)
    {ind.tri<-ifelse(all(c(alpha,beta,gamma)>=0), TRUE, FALSE)
    }  else
    {ind.tri<-ifelse(all(c(alpha,beta,gamma)>0), TRUE, FALSE)}
  }
  list(
   in.tri=ind.tri,
    barycentric=c(alpha,beta,gamma)
  )
} #end of the function
#'

#################################################################

#' @title Check whether all points in a data set are inside the triangle
#'
#' @description Checks if all the data points in the 2D data set, \code{Xp}, lie in the triangle, \code{tri},
#' using the barycentric coordinates, generally denoted as \eqn{(\alpha,\beta,\gamma)}.
#'
#' If all (normalized or non-normalized) barycentric coordinates of a point are positive then the point is
#' inside the triangle, if all are nonnegative with one or more are zero, then the point falls in the boundary.
#' If some of the barycentric coordinates are negative, then the point falls outside the triangle.
#'
#' \code{boundary} is a logical argument (default=\code{TRUE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if a point lies in the closure of the triangle (i.e., interior and boundary combined);
#' else, it checks if the point lies in the interior of the triangle.
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param boundary A logical parameter (default=\code{FALSE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if a point lies in the closure of the triangle (i.e., interior and boundary combined)
#' else it checks if the point lies in the interior of the triangle.
#'
#' @return A logical output, if all data points in \code{Xp} are inside the triangle, \code{tri}, the output is \code{TRUE},
#' else it is \code{FALSE}.
#'
#' @seealso \code{\link{in.triangle}} and \code{\link[interp]{on.convex.hull}} from the \code{interp} package for documentation for \code{in.convex.hull}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2); p<-c(1.4,1.2)
#'
#' Tr<-rbind(A,B,C)
#'
#' in.tri.all(p,Tr)
#'
#' #for the vertex A
#' in.tri.all(A,Tr)
#' in.tri.all(A,Tr,boundary=FALSE)
#'
#' #for a point on the edge AB
#' D3<-(A+B)/2
#' in.tri.all(D3,Tr)
#' in.tri.all(D3,Tr,boundary=FALSE)
#'
#' #data set
#' n<-10
#' Xp<-cbind(runif(n),runif(n))
#' in.tri.all(Xp,Tr,boundary=TRUE)
#'
#' Xp<-runif.std.tri(n)$gen.points
#' in.tri.all(Xp,Tr)
#' in.tri.all(Xp,Tr,boundary=FALSE)
#'
#' Xp<-runif.tri(n,Tr)$g
#' in.tri.all(Xp,Tr)
#' in.tri.all(Xp,Tr,boundary=FALSE)
#' }
#'
#' @export in.tri.all
in.tri.all <- function(Xp,tri,boundary=TRUE)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { intri<-in.triangle(Xp,tri,boundary)$in.tri
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  n<-nrow(Xp)
  cnt<-0; i<-1; intri<-TRUE
  while (i<=n & cnt==0)
  { pnt<-Xp[i,]
  if (!in.triangle(pnt,tri,boundary)$in.tri)
  {cnt<-1; intri<-FALSE}
  else
  {i<-i+1}
  }
  }
 intri
} #end of the function
#'

#################################################################

#' @title Check whether a triangle is a standard equilateral triangle
#'
#' @description Checks whether the triangle, \code{tri}, is the standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' or not.
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return \code{TRUE} if \code{tri} is a standard equilateral triangle, else \code{FALSE}.
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)  #try adding +10^(-16) to each vertex
#' is.std.eq.tri(Te)
#'
#' is.std.eq.tri(rbind(B,C,A))
#'
#' Tr<-rbind(A,B,-C)
#' is.std.eq.tri(Tr)
#'
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' is.std.eq.tri(Tr)
#' }
#'
#' @export is.std.eq.tri
is.std.eq.tri <- function(tri)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('the argument must be numeric and of dimension 3x2')}

  p1<-as.numeric(tri[1,]); p2<-as.numeric(tri[2,]); p3<-as.numeric(tri[3,]);
  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2)

  checkA<-isTRUE(all.equal(p1,A)) + isTRUE(all.equal(p2,A)) + isTRUE(all.equal(p3,A))
  checkB<-isTRUE(all.equal(p1,B)) + isTRUE(all.equal(p2,B)) + isTRUE(all.equal(p3,B))
  checkC<-isTRUE(all.equal(p1,C)) + isTRUE(all.equal(p2,C)) + isTRUE(all.equal(p3,C))

  eq.tri<-sum(checkA+checkB+checkC)==3 & min(checkA,checkB,checkC)==1
  eq.tri
} #end of the function
#'


#################################################################

#' @title The index of the \eqn{CM}-vertex region in a triangle that contains a given point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' the triangle \code{tri}\eqn{=(y_1,y_2,y_3)} with vertex regions are constructed with center of mass \eqn{CM=(y_1+y_2+y_3)/3}
#' (see the plots in the example for illustrations).
#'
#' The vertices of triangle, \code{tri}, are labeled as \eqn{1,2,3}
#' according to the row number the vertex is recorded in \code{tri}. If the point, \code{p}, is not inside \code{tri}, then the
#' function yields \code{NA} as output. The corresponding vertex region is the polygon with the vertex, \eqn{CM}, and
#' midpoints of the edges adjacent to the vertex.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds})
#'
#' @param p A 2D point for which \eqn{CM}-vertex region it resides in is to be determined in the triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the \eqn{CM}-vertex region that contains point, \code{p} in the triangle \code{tri}.}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rv.tri.cent}}, \code{\link{rv.triCC}}, \code{\link{rv.bas.triCM}},
#' \code{\link{rv.bas.triCC}}, \code{\link{rv.bas.tri.cent}}, and \code{\link{rvTeCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.6,2);
#' Tr<-rbind(A,B,C);
#' P<-c(1.4,1.2)
#' rv.triCM(P,Tr)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.tri(n,Tr)$g
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rv.triCM(Xp[i,],Tr)$rv)
#' Rv
#'
#' CM<-(A+B+C)/3
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".")
#' L<-Ds; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(Rv))
#'
#' txt<-rbind(Tr,CM,D1,D2,D3)
#' xc<-txt[,1]+c(-.02,.02,.02,-.02,.02,-.01,-.01)
#' yc<-txt[,2]+c(-.02,-.04,.06,-.02,.02,.06,-.06)
#' txt.str<-c("rv=1","rv=2","rv=3","CM","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export rv.triCM
rv.triCM <- function(p,tri)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) )
  {stop('tri must be numeric')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (in.triangle(p,tri,boundary=TRUE)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];
    a1<-y1[1]; a2<-y1[2]; b1<-y2[1]; b2<-y2[2]; c1<-y3[1]; c2<-y3[2];

    CM<-(y1+y2+y3)/3

    D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
    x<-p[1]; y<-p[2];

    if (in.triangle(p,rbind(y1,D3,CM),boundary=TRUE)$in.tri | in.triangle(p,rbind(y1,CM,D2),boundary=TRUE)$in.tri)
    {rv<-1}
    else
    {
      if (in.triangle(p,rbind(D3,y2,CM),boundary=TRUE)$in.tri | in.triangle(p,rbind(y2,D1,CM),boundary=TRUE)$in.tri)
      {rv<-2}
      else
      {rv<-3}
    }
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The labels of the vertices of a triangle in the basic triangle form
#'
#' @description Labels the vertices of triangle, \code{tri}, as \eqn{ABC} so that  \eqn{AB} is the longest edge,  \eqn{BC} is the second longest
#' and  \eqn{AC} is the shortest edge (the order is as in the basic triangle).
#'
#' The standard basic triangle form is \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} where \eqn{c_1} is in \eqn{[0,1/2]},
#' \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' The option \code{scaled} a logical argument for scaling the resulting triangle or not.
#' If \code{scaled=TRUE}, then the resulting triangle is scaled to be a regular basic triangle, i.e., longest edge having unit length,
#' else (i.e., if \code{scaled=FALSE} which is the default), the new triangle \eqn{T(A,B,C)} is nonscaled,
#' i.e., the longest edge \eqn{AB} may not be of unit length.
#' The vertices of the resulting triangle (whether scaled or not) is presented in the order of vertices of the corresponding
#' basic triangle, however when scaled the triangle is equivalent to the basic triangle \eqn{T_b} up to translation and rotation.
#' That is, this function converts any triangle to a basic triangle (up to translation and rotation),
#' so that the output triangle is $T(A',B',C')$ so that edges in decreasing length are $A'B'$, $B'C'$, and $A'C'$.
#' Most of the times,
#' the resulting triangle will still need to be translated and/or rotated to be in the standard basic triangle form.
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param scaled A logical argument for scaling the resulting basic triangle. If \code{scaled=TRUE}, then the resulting triangle is
#' scaled to be a regular basic triangle, i.e., longest edge having unit length,
#' else the new triangle \eqn{T(A,B,C)} is nonscaled. The default is \code{scaled=FALSE}.
#'
#' @return A \code{list} with three elements
#' \item{tri}{The vertices of the basic triangle stacked row-wise and labeled row-wise as \eqn{A}, \eqn{B}, \eqn{C}.}
#' \item{desc}{Description of the edges based on the vertices, i.e.,
#'  \code{"Edges (in decreasing length are) AB, BC, and AC"}.}
#' \item{orig.order}{Row order of the input triangle, \code{tri},
#' when converted to the scaled version of the basic triangle}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#'
#' as.bas.tri(rbind(A,B,C))
#' as.bas.tri(rbind(B,C,A))
#'
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' as.bas.tri(rbind(A,B,C))
#' as.bas.tri(rbind(A,C,B))
#' as.bas.tri(rbind(B,A,C))
#' }
#'
#' @export
as.bas.tri <- function(tri,scaled=FALSE)
{
  tri<-as.matrix(tri)#
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('the argument must be numeric and of dimension 3x2')}

  v1<-tri[1,];v2<-tri[2,];v3<-tri[3,];
  de1<-Dist(v1,v2); de2<-Dist(v2,v3); de3<-Dist(v1,v3);
  max.de =max(de1,de2,de3) #maximum edge length of the triangle
  ord<-order(c(de1,de2,de3),decreasing = TRUE)
  vord<-cbind(c(1,2,1),c(2,3,3))
  vord<-vord[ord,]
  C<-tri[setdiff(1:3,vord[1,]),]
  A<-tri[setdiff(1:3,vord[2,]),]
  B<-tri[setdiff(1:3,vord[3,]),]

  orig.ord<-c(setdiff(1:3,vord[2,]),setdiff(1:3,vord[3,]),setdiff(1:3,vord[1,]))
  ifelse(scaled, Tr<-rbind(A,B,C)/max.de, Tr<-rbind(A,B,C))
  row.names(Tr)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edges (in decreasing length are) AB, BC, and AC"

  list(tri=Tr, #vertex labeling
       desc=edge.desc,
       orig.order=orig.ord #order of vertices in argument tri
  )
} #end of the function
#'

#################################################################

#' @title Converting a triangle to the standard basic triangle form form
#'
#' @description This function transforms any triangle, \code{tri}, to the standard basic triangle form.
#'
#' The standard basic triangle form is \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} where \eqn{c_1} is in \eqn{[0,1/2]},
#' \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle form by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with two elements
#' \item{Cvec}{The nontrivial vertex \eqn{C=(c_1,c_2)} in the standard basic triangle form \eqn{T_b}.}
#' \item{orig.order}{Row order of the input triangle, \code{tri},
#' when converted to the standard basic triangle form \eqn{T_b}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#'
#' tri2std.bas.tri(rbind(A,B,C))
#' tri2std.bas.tri(rbind(B,C,A))
#'
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' tri2std.bas.tri(rbind(A,B,C))
#' tri2std.bas.tri(rbind(A,C,B))
#' tri2std.bas.tri(rbind(B,A,C))
#' }
#'
#' @export
tri2std.bas.tri <- function(tri)
{
  tri.scaled <-as.bas.tri(tri,scaled = TRUE) # scaling so that the longest edge is of unit length
  tri.new =tri.scaled$tri

  A<-tri.new[1,];B<-tri.new[2,];C<-tri.new[3,];
  ac<-Dist(A,C); bc<-Dist(B,C);

  ab.vec = B-A
  ac.vec = C-A
  bc.vec = C-B

  theta1 <- acos( sum(ab.vec*ac.vec) / ( sqrt(sum(ab.vec * ab.vec)) * sqrt(sum(ac.vec * ac.vec)) ) ) #angle between edges AB and AC
  theta2 <- acos( sum(ab.vec*bc.vec) / ( sqrt(sum(ab.vec * ab.vec)) * sqrt(sum(bc.vec * bc.vec)) ) ) #angle between edges AB and BC

  c1=ac*cos(theta1)
  c2=bc*sin(theta2)

  list(Cvec=c(c1,c2), #the nontrivial vertex C in the standard basic triangle form
       orig.order=tri.scaled$orig.ord #order of vertices in argument tri
  )
} #end of the function
#'

#################################################################

# funsCartBary
#'
#' @title Converts of a point in Cartesian coordinates to Barycentric coordinates and vice versa
#'
#' @description
#' Two functions: \code{cart2bary} and \code{bary2cart}.
#'
#' \code{cart2bary} converts Cartesian coordinates of a given point \code{P}\eqn{=(x,y)} to barycentric coordinates
#' (in the normalized form) with respect to the triangle \code{tri}\eqn{=(v_1,v_2,v_3)}
#' with vertex labeling done row-wise in \code{tri} (i.e., row \eqn{i} corresponds to vertex \eqn{v_i} for \eqn{i=1,2,3}).
#'
#' \code{bary2cart} converts barycentric coordinates of the point \code{P}\eqn{=(t_1,t_2,t_3)} (not necessarily normalized) to
#' Cartesian coordinates according to the coordinates of the triangle, \code{tri}.
#' For information on barycentric coordinates, see (\insertCite{weisstein-barycentric;textual}{pcds}).
#'
#' @param P A 2D point for \code{cart2bary}, and a \code{vector} of three \code{numeric} entries for \code{bary2cart}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return \code{cart2bary} returns the barycentric coordinates of a given point \code{P}\eqn{=(x,y)} and
#' \code{bary2cart} returns the Cartesian coordinates of the point \code{P}\eqn{=(t_1,t_2,t_3)} (not necessarily normalized)
#'
#' @name funsCartBary
NULL
#'
#' @rdname funsCartBary
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for cart2bary
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tr<-rbind(A,B,C)
#'
#' cart2bary(A,Tr)
#' cart2bary(c(.3,.2),Tr)
#' }
#'
#' @export cart2bary
cart2bary <- function(P,tri)
{
  if (!is.point(P))
  {stop('P must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  v1<-tri[1,];v2<-tri[2,];v3<-tri[3,];
  T<-cbind(v1,v2)-cbind(v3,v3)
  lam<-solve(T) %*% (matrix(P-v3))
  c(lam,1-sum(lam))
} #end of the function
#'
#' @rdname funsCartBary
#'
#' @examples
#' \dontrun{
#' #Examples for bary2cart
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tr<-rbind(A,B,C)
#'
#' bary2cart(c(.3,.2,.5),Tr)
#' bary2cart(c(6,2,4),Tr)
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export bary2cart
bary2cart <- function(P,tri)
{
  if (!is.point(P,3))
  {stop('P must be a numeric 3D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  P<-P/sum(P)  #normalized barycentric coordinates
  R<-t(tri)
  cc<-R %*% matrix(P)  #Cartesian coordinates
  as.vector(cc)
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CM}-vertex region in a standard basic triangle form that contains a point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' the standard basic triangle form \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0}
#' and \eqn{(1-c_1)^2+c_2^2 \le 1} and vertex regions are based on the center of mass CM=((1+c1)/3,c2/3) of \eqn{T_b}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the standard basic triangle form \eqn{T_b} are labeled as
#' \eqn{1=(0,0)}, \eqn{2=(1,0)},and \eqn{3=(c_1,c_2)} also according to the row number the vertex is recorded in \eqn{T_b}.
#' If the point, \code{p}, is not inside \eqn{T_b}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon with the vertex, \eqn{CM}, and
#' midpoints of the edges adjacent to the vertex.
#'
#' Any given triangle can be mapped to the standard basic triangle form
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-PE;textual}{pcds})
#'
#' @param p A 2D point for which \eqn{CM}-vertex region it resides in is to be determined in the
#' standard basic triangle form \eqn{T_b}.
#' @param c1,c2 Positive real numbers which constitute the upper vertex of the standard basic triangle form (i.e., the vertex
#' adjacent to the shorter edges of \eqn{T_b}); \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and
#' \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the \eqn{CM}-vertex region that contains point, \code{p} in the standard basic triangle form \eqn{T_b}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}
#' with row \eqn{1=(0,0)}, row \eqn{2=(1,0)}, and row \eqn{3=(c_1,c_2)}.}
#'
#' @seealso \code{\link{rv.triCM}}, \code{\link{rv.tri.cent}}, \code{\link{rv.triCC}},
#' \code{\link{rv.bas.triCC}}, \code{\link{rv.bas.tri.cent}}, and \code{\link{rvTeCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' #' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6
#' P<-c(.4,.2)
#' rv.bas.triCM(P,c1,c2)
#'
#' A<-c(0,0);B<-c(1,0);C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#' CM<-(A+B+C)/3
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rv.bas.triCM(Xp[i,],c1,c2)$rv)
#' Rv
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,xlab="",ylab="",axes="T",pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(Xp,pch=".")
#' polygon(Tb)
#' L<-Ds; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(Rv))
#'
#' txt<-rbind(Tb,CM,Ds)
#' xc<-txt[,1]+c(-.03,.03,.02,-.01,.06,-.05,.0)
#' yc<-txt[,2]+c(.02,.02,.02,.04,.02,.02,-.03)
#' txt.str<-c("A","B","C","CM","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' plot(Tb,xlab="",ylab="",axes="T",pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-Ds; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' RV1<-(A+D3+CM+D2)/4
#' RV2<-(B+D3+CM+D1)/4
#' RV3<-(C+D2+CM+D1)/4
#'
#' txt<-rbind(RV1,RV2,RV3)
#' xc<-txt[,1]
#' yc<-txt[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(Tb,CM,Ds)
#' xc<-txt[,1]+c(-.03,.03,.02,-.01,.04,-.03,.0)
#' yc<-txt[,2]+c(.02,.02,.02,.04,.02,.02,-.03)
#' txt.str<-c("A","B","C","CM","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export rv.bas.triCM
rv.bas.triCM <- function(p,c1,c2)
{
  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  p1<-c(0,0); p2<-c(1,0); p3<-c(c1,c2);
  Tb<-rbind(p1,p2,p3)

  if (in.triangle(p,Tb,boundary = TRUE)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    x<-p[1]; y<-p[2];

    if (c1<1/2)
    {
      if ( y<=c2*x/(1+c1) && y>=c2 * (-1 + 2 * x) / (2 * c1 - 1))
      {rv<-2}
      else
      {
        if (y>c2*x/(1+c1) && y>=c2 * (-1 + x) / (c1 - 2))
        {rv<-3}
        else
        {rv<-1}
      }
    }
    else
    {
      if (c1 == 1/2)
      {
        if ( y<=c2*x/(1+c1) && x>=c1)
        {rv<-2}
        else
        {
          if (y>c2*x/(1+c1) && y>=c2 * (-1 + x) / (c1 - 2))
          {rv<-3}
          else
          {rv<-1}
        }
      }
      else
      {if ( y<=c2*x/(1+c1) && y<=c2 * (-1 + 2 * x) / (2 * c1 - 1))
      {rv<-2}
        else
        {
          if (y>c2*x/(1+c1) && y>=c2 * (-1 + x) / (c1 - 2))
          {rv<-3}
          else
          {rv<-1}
        }
      }
    }
  }
  row.names(Tb)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=Tb #vertex labeling
  )
} #end of the function
#'



#################################################################

#' @title The index of the \eqn{CM}-edge region in a triangle that contains the point
#'
#' @description Returns the index of the edge whose region contains point, \code{p}, in
#' the triangle \code{tri}\eqn{=T(A,B,C)} with edge regions based on center of mass \eqn{CM=(A+B+C)/3}.
#'
#' Edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' Edge region 1 is the triangle \eqn{T(B,C,CM)}, edge region 2 is \eqn{T(A,C,CM)}, and
#' edge region 3 is \eqn{T(A,B,CM)}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CM}-edge region it resides in is to be determined in the triangle
#' \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with three elements
#' \item{re}{Index of the \eqn{CM}-edge region that contains point, \code{p} in the triangle \code{tri}.}
#' \item{tri}{The vertices of the triangle, where row labels are \eqn{A}, \eqn{B}, and \eqn{C}
#' with edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.}
#' \item{desc}{Description of the edge labels}
#'
#' @seealso \code{\link{re.tri.cent}}, \code{\link{re.bas.triCM}}, \code{\link{re.bas.tri.cent}},
#' \code{\link{reTeCM}}, and \code{\link{edge.reg.triCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' P<-c(1.4,1.2)
#' re.triCM(P,Tr)
#'
#' P<-c(1.5,1.61)
#' re.triCM(P,Tr)
#'
#' CM<-(A+B+C)/3
#'
#' n<-20  #try also n<-40
#' Xp<-runif.tri(n,Tr)$g
#'
#' re<-vector()
#' for (i in 1:n)
#'   re<-c(re,re.triCM(Xp[i,],Tr)$re)
#' re
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(Xp,pch=".")
#' polygon(Tr)
#' L<-Tr; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(re))
#'
#' txt<-rbind(Tr,CM)
#' xc<-txt[,1]
#' yc<-txt[,2]
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#'
#' p1<-(A+B+CM)/3
#' p2<-(B+C+CM)/3
#' p3<-(A+C+CM)/3
#'
#' plot(Tr,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-Tr; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tr,CM,p1,p2,p3)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.02,.02,.02)
#' yc<-txt[,2]+c(.02,.02,.04,.05,.02,.02,.02)
#' txt.str<-c("A","B","C","CM","re=3","re=1","re=2")
#' text(xc,yc,txt.str)
#' }
#'
#' @export re.triCM
re.triCM <- function(p,tri)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (in.triangle(p,tri,boundary=TRUE)$in.tri==FALSE)
  {reled<-NA
  } else
  {
    y1<-tri[1,];y2<-tri[2,]; y3<-tri[3,];
    CM<-(1/3)*(y1+y2+y3);

    if (in.triangle(p,rbind(y1,y2,CM),boundary = TRUE)$in.tri)
    {reled<-3
    } else if (in.triangle(p,rbind(y2,y3,CM),boundary = TRUE)$in.tri)
    {reled<-1
    } else
    {reled<-2}
  }
  row.names(tri)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=reled, #relative edge
       tri=tri, #vertex labeling
       desc=edge.desc
  )
} #end of the function
#'

#################################################################

#' @title The vertices of the \eqn{CM}-edge region in a triangle that contains the point
#'
#' @description Returns the edge whose region contains point, \code{p}, in
#' the triangle \code{tri}\eqn{=T(A,B,C)} with edge regions based on center of mass \eqn{CM=(A+B+C)/3}.
#'
#' This function is related to \code{\link{re.triCM}}, but unlike \code{\link{re.triCM}} the related edges are given as vertices \code{ABC}
#' for \eqn{re=3}, as \code{BCA} for \eqn{re=1} and as \code{CAB} for \eqn{re=2} where edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC},
#' and 2 for edge \eqn{AC}. The vertices are given one vertex in each row in the output, e.g., \eqn{ABC} is printed
#' as \code{rbind(A,B,C)}, where row 1 has the entries of vertex A, row 2 has the entries of vertex B,
#' and row 3 has the entries of vertex C.
#'
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#'
#' Edge region for BCA is the triangle \eqn{T(B,C,CM)}, edge region CAB is \eqn{T(A,C,CM)}, and edge region ABC is \eqn{T(A,B,CM)}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CM}-edge region it resides in is to be determined in the triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return The \eqn{CM}-edge region that contains point, \code{p} in the triangle \code{tri}. The related edges are given as
#' vertices \code{ABC} for \eqn{re=3}, as \code{BCA} for \eqn{re=1} and as \code{CAB} for \eqn{re=2} where edges are labeled as 3
#' for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.
#'
#' @seealso \code{\link{re.tri.cent}}, \code{\link{re.triCM}}, \code{\link{re.bas.triCM}},
#' \code{\link{re.bas.tri.cent}}, \code{\link{reTeCM}}, and \code{\link{edge.reg.triCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' P<-c(.4,.2)  #try also P<-as.numeric(runif.tri(1,Tr)$g)
#' edge.reg.triCM(P,Tr)
#'
#' P<-c(1.8,.5)
#' edge.reg.triCM(P,Tr)
#'
#' CM<-(A+B+C)/3
#' p1<-(A+B+CM)/3
#' p2<-(B+C+CM)/3
#' p3<-(A+C+CM)/3
#'
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-Tr; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tr,CM,p1,p2,p3)
#' xc<-txt[,1]+c(-.02,.02,.02,-.05,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,.02,0,0,0)
#' txt.str<-c("A","B","C","CM","re=T(A,B,CM)","re=T(B,C,CM)","re=T(A,C,CM)")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
edge.reg.triCM <- function(p,tri)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];

  if (in.triangle(p,tri)$in.tri==FALSE)
  {reled<-NA;return(reled);stop('point is not inside the triangle')}

  CM<-(y1+y2+y3)/3;
  M3<-(y1+y2)/2; M2<-(y1+y3)/2; M1<-(y2+y3)/2;
  x<-p[1]; y<-p[2];
  ifelse(as.numeric(abs(Line(y1,CM,x)$sl))==Inf,cond31<-sign(x-y1[1])==sign(M3[1]-y1[1]),
         cond31<-sign(y-Line(y1,CM,x)$y)==sign(M3[2]-Line(y1,CM,M3[1])$y))

  ifelse(as.numeric(abs(Line(y2,CM,x)$sl))==Inf,cond32<-sign(x-y2[1])==sign(M3[1]-y2[1]),
         cond32<-sign(y-Line(y2,CM,x)$y)==sign(M3[2]-Line(y2,CM,M3[1])$y))

  ifelse(as.numeric(abs(Line(y2,CM,x)$sl))==Inf,cond11<-sign(x-y2[1])==sign(M1[1]-y2[1]),
         cond11<-sign(y-Line(y2,CM,x)$y)==sign(M1[2]-Line(y2,CM,M1[1])$y))

  ifelse(as.numeric(abs(Line(y3,CM,x)$sl))==Inf,cond13<-sign(x-y3[1])==sign(M1[1]-y3[1]),
         cond13<-sign(y-Line(y3,CM,x)$y)==sign(M1[2]-Line(y3,CM,M1[1])$y))

  if (cond31 && cond32)
    reled<-rbind(y1,y2,CM)
  else
  {
    if ( cond11 && cond13)
      reled<-rbind(y2,y3,CM)
    else reled<-rbind(y3,y1,CM)
  }
  reled
} #end of the function
#'

#################################################################

#' @title The index of the edge region in a triangle that contains the point
#'
#' @description Returns the index of the edge whose region contains point, \code{p}, in
#' the triangle \code{tri}\eqn{=T(A,B,C)} with edge regions based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}.
#'
#' Edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' Edge region 1 is the triangle \eqn{T(B,C,M)}, edge region 2 is \eqn{T(A,C,M)}, and edge region 3 is \eqn{T(A,B,M)}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p A 2D point for which \code{M}-edge region it resides in is to be determined in the triangle
#' \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}.
#'
#' @return A \code{list} with three elements
#' \item{re}{Index of the \code{M}-edge region that contains point, \code{p} in the triangle \code{tri}.}
#' \item{tri}{The vertices of the triangle, where row labels are \eqn{A}, \eqn{B}, and \eqn{C}
#' with edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.}
#' \item{desc}{Description of the edge labels}
#'
#' @seealso \code{\link{re.triCM}}, \code{\link{re.bas.triCM}},
#' \code{\link{re.bas.tri.cent}}, \code{\link{reTeCM}}, and \code{\link{edge.reg.triCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' P<-c(1.4,1.2)
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' re.tri.cent(P,Tr,M)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.tri(n,Tr)$g
#'
#' re<-vector()
#' for (i in 1:n)
#'   re<-c(re,re.tri.cent(Xp[i,],Tr,M)$re)
#' re
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#'
#' plot(Tr,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".")
#' L<-Tr; R<-rbind(M,M,M)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(re))
#'
#' txt<-rbind(Tr,M)
#' xc<-txt[,1]
#' yc<-txt[,2]
#' txt.str<-c("A","B","C","M")
#' text(xc,yc,txt.str)
#'
#' p1<-(A+B+M)/3
#' p2<-(B+C+M)/3
#' p3<-(A+C+M)/3
#'
#' plot(Tr,xlab="",ylab="", main="Illustration of M-edge regions in a triangle",
#' axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-Tr; R<-rbind(M,M,M)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tr,M,p1,p2,p3)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.02,.02,.02)
#' yc<-txt[,2]+c(.02,.02,.04,.05,.02,.02,.02)
#' txt.str<-c("A","B","C","M","re=3","re=1","re=2")
#' text(xc,yc,txt.str)
#' }
#'
#' @export re.tri.cent
re.tri.cent <- function(p,tri,M)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {
    M<-bary2cart(M,tri)
  }
  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (in.triangle(p,tri,boundary=TRUE)$in.tri==FALSE)
  {reled<-NA
  } else
  {
    y1<-tri[1,];y2<-tri[2,]; y3<-tri[3,];

    if (in.triangle(p,rbind(y1,y2,M),boundary = TRUE)$in.tri)
      reled<-3
    else
    {
      if (in.triangle(p,rbind(y2,y3,M),boundary = TRUE)$in.tri)
        reled<-1
      else reled<-2
    }
  }
  row.names(tri)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=reled, #relative edge
       tri=tri, #vertex labeling
       desc=edge.desc
  )
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CM}-edge region in a standard basic triangle form that contains a point
#'
#' @description Returns the index of the edge whose region contains point, \code{p}, in the
#' standard basic triangle form \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2)} where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1} with
#' edge regions based on center of mass \eqn{CM=(A+B+C)/3}.
#'
#' Edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' Edge region 1 is the triangle \eqn{T(B,C,CM)}, edge region 2 is \eqn{T(A,C,CM)}, and edge region 3 is \eqn{T(A,B,CM)}.
#'
#' Any given triangle can be mapped to the standard basic triangle form
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CM}-edge region it resides in is to be determined in the
#' standard basic triangle form \eqn{T_b}.
#' @param c1,c2 Positive real numbers which constitute the upper vertex of the standard basic triangle form (i.e., the vertex
#' adjacent to the shorter edges of \eqn{T_b}); \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and
#' \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' @return A \code{list} with three elements
#'
#' \item{re}{Index of the \eqn{CM}-edge region that contains point, \code{p} in the standard basic triangle form \eqn{T_b}}
#' \item{tri}{The vertices of the triangle, where row labels are \eqn{A=(0,0)}, \eqn{B=(1,0)}, and \eqn{C=(c_1,c_2)}
#' with edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.}
#' \item{desc}{Description of the edge labels}
#'
#' @seealso \code{\link{re.triCM}}, \code{\link{re.tri.cent}},
#' \code{\link{re.bas.tri.cent}}, \code{\link{reTeCM}}, and \code{\link{edge.reg.triCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6
#' P<-c(.4,.2)
#' re.bas.triCM(P,c1,c2)
#'
#' A<-c(0,0);B<-c(1,0);C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#' CM<-(A+B+C)/3
#'
#' re.bas.triCM(A,c1,c2)
#' re.bas.triCM(B,c1,c2)
#' re.bas.triCM(C,c1,c2)
#' re.bas.triCM(CM,c1,c2)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' re<-vector()
#' for (i in 1:n)
#'   re<-c(re,re.bas.triCM(Xp[i,],c1,c2)$re)
#' re
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(Xp,pch=".")
#' polygon(Tb)
#' L<-Tb; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(re))
#'
#' txt<-rbind(Tb,CM)
#' xc<-txt[,1]+c(-.03,.03,.02,0)
#' yc<-txt[,2]+c(.02,.02,.02,-.04)
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#' }
#'
#' @export re.bas.triCM
re.bas.triCM <- function(p,c1,c2)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2);
  Tb<-rbind(y1,y2,y3)

  if (in.triangle(p,Tb,boundary = TRUE)$in.tri==FALSE)
  {reled<-NA
  } else
  {
    CM<-(1/3)*(y1+y2+y3);

    if (in.triangle(p,rbind(y1,y2,CM),boundary = TRUE)$in.tri)
      reled<-3
    else
    {
      if (in.triangle(p,rbind(y2,y3,CM),boundary = TRUE)$in.tri)
        reled<-1
      else reled<-2
    }
  }
  row.names(Tb)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=reled, #relative edge
       tri=Tb, #vertex labeling
       desc=edge.desc #description of the edge labels
  )
} #end of the function
#'

#################################################################

#' @title The index of the edge region in a standard basic triangle form that contains a point
#'
#' @description Returns the index of the edge whose region contains point, \code{p}, in
#' the standard basic triangle form \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2))} and edge regions based on center
#' \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the standard basic triangle form \eqn{T_b}.
#'
#' Edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' Edge region 1 is the triangle \eqn{T(B,C,M)}, edge region 2 is \eqn{T(A,C,M)}, and edge region 3 is \eqn{T(A,B,M)}.
#' In the standard basic triangle form \eqn{T_b} \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle form by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p A 2D point for which \code{M}-edge region it resides in is to be determined in the
#' standard basic triangle form \eqn{T_b}.
#' @param c1,c2 Positive real numbers which constitute the upper vertex of the standard basic triangle form (i.e., the vertex
#' adjacent to the shorter edges of \eqn{T_b}); \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and
#' \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle form \eqn{T_b}.
#'
#' @return A \code{list} with three elements
#'
#' \item{re}{Index of the \code{M}-edge region that contains point, \code{p} in the standard basic triangle form \eqn{T_b}.}
#' \item{tri}{The vertices of the triangle, where row labels are \eqn{A}, \eqn{B}, and \eqn{C}
#' with edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.}
#' \item{desc}{Description of the edge labels}
#'
#' @seealso \code{\link{re.triCM}}, \code{\link{re.tri.cent}},
#' \code{\link{re.bas.tri.cent}}, \code{\link{reTeCM}}, and \code{\link{edge.reg.triCM}}
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
#' M<-c(.6,.2)
#'
#' P<-c(.4,.2)
#' re.bas.tri.cent(P,c1,c2,M)
#'
#' A<-c(0,0);B<-c(1,0);C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' re<-vector()
#' for (i in 1:n)
#'   re<-c(re,re.bas.tri.cent(Xp[i,],c1,c2,M)$re)
#' re
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(Xp,pch=".")
#' polygon(Tb)
#' L<-Tb; R<-rbind(M,M,M)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(re))
#'
#' txt<-rbind(Tb,M)
#' xc<-txt[,1]+c(-.03,.03,.02,0)
#' yc<-txt[,2]+c(.02,.02,.02,-.03)
#' txt.str<-c("A","B","C","M")
#' text(xc,yc,txt.str)
#' }
#'
#' @export re.bas.tri.cent
re.bas.tri.cent <- function(p,c1,c2,M)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 arguments must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2);
  Tb<-rbind(y1,y2,y3)

  if (dimension(M)==3)
  {
    M<-bary2cart(M,Tb)
  }

  if (in.triangle(M,Tb,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (in.triangle(p,Tb,boundary = TRUE)$in.tri==FALSE)
  {reled<-NA
  } else
  {

    if (in.triangle(p,rbind(y1,y2,M),boundary = TRUE)$in.tri)
      reled<-3
    else
    {
      if (in.triangle(p,rbind(y2,y3,M),boundary = TRUE)$in.tri)
        reled<-1
      else reled<-2
    }
  }
  row.names(Tb)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=reled, #relative edge
       tri=Tb, #vertex labeling
       desc=edge.desc
  )
} #end of the function
#'

#################################################################

#' @title The index of the edge region in the standard equilateral triangle that contains a point
#'
#' @description Returns the index of the edge whose region contains point, \code{p}, in
#' the standard equilateral triangle \eqn{T_e=T(A=(0,0),B=(1,0),C=(1/2,\sqrt{3}/2))} with edge regions based on
#' center of mass \eqn{CM=(A+B+C)/3}.
#'
#' Edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' Edge region 1 is the triangle \eqn{T(B,C,M)}, edge region 2 is \eqn{T(A,C,M)}, and edge region 3 is \eqn{T(A,B,M)}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CM}-edge region it resides in is to be determined in the
#' the standard equilateral triangle \eqn{T_e}.
#'
#' @return A \code{list} with three elements
#' \item{re}{Index of the \eqn{CM}-edge region that contains point, \code{p} in the standard equilateral triangle \eqn{T_e}}
#' \item{tri}{The vertices of the standard equilateral triangle \eqn{T_e}, where row labels are \eqn{A}, \eqn{B}, and \eqn{C}
#' with edges are labeled as 3 for edge \eqn{AB}, 1 for edge \eqn{BC}, and 2 for edge \eqn{AC}.}
#' \item{desc}{Description of the edge labels}
#'
#' @seealso \code{\link{re.triCM}}, \code{\link{re.tri.cent}},
#' \code{\link{re.bas.triCM}}, \code{\link{re.bas.tri.cent}}, and \code{\link{edge.reg.triCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P<-c(.4,.2)
#' reTeCM(P)
#'
#' A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' CM<-(A+B+C)/3
#'
#' n<-20  #try also n<-40
#' Xp<-runif.std.tri(n)$gen.points
#'
#' re<-vector()
#' for (i in 1:n)
#'   re<-c(re,reTeCM(Xp[i,])$re)
#' re
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,asp=1,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' points(Xp,pch=".")
#' polygon(Te)
#' L<-Te; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(re))
#'
#' txt<-rbind(Te,CM)
#' xc<-txt[,1]+c(-.03,.03,.03,-.06)
#' yc<-txt[,2]+c(.02,.02,.02,.03)
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#'
#' p1<-(A+B+CM)/3
#' p2<-(B+C+CM)/3
#' p3<-(A+C+CM)/3
#'
#' plot(Te,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' L<-Te; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Te,CM,p1,p2,p3)
#' xc<-txt[,1]+c(-.03,.03,.03,-.06,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,.03,0,0,0)
#' txt.str<-c("A","B","C","CM","re=3","re=1","re=2")
#' text(xc,yc,txt.str)
#' }
#'
#' @export reTeCM
reTeCM <- function(p)
{
  if (!is.point(p))
  {stop('the argument must be a numeric 2D point')}

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
  Te<-rbind(A,B,C)
  if (in.triangle(p,Te,boundary=TRUE)$in.tri==FALSE)
  {reled<-NA
  } else
  {
    if (p[2]<= .5773502693*p[1] && p[2]<= .5773502693-.5773502693*p[1])
      reled<-3
    else
    {
      reled<-ifelse (p[1]>=1/2, 1, 2)
    }
  }
  row.names(Te)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=reled, #relative edge
       tri=Te, #vertex labeling
       desc=edge.desc
  )
} #end of the function
#'

#################################################################

#' @title The indices of the \eqn{CM}-edge regions in a triangle that contains the points in a give data set
#'
#' @description Returns the indices of the edges whose regions contain the points in data set \code{Xp} in
#' a triangle \code{tri}\eqn{=(A,B,C)} and edge regions are based on the center of mass \eqn{CM} of \code{tri}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the  triangle \code{tri} are labeled as
#' \eqn{1=A}, \eqn{2=B}, and \eqn{3=C} also according to the row number the vertex is recorded in \code{tri}
#' and the corresponding edges are \eqn{1=BC}, \eqn{2=AC}, and \eqn{3=AB}.
#'
#' If a point in \code{Xp} is not inside \code{tri}, then the function yields \code{NA} as output for that entry.
#' The corresponding edge region is the polygon
#' with the vertex, \eqn{CM}, and vertices other than the non-adjacent vertex, i.e., edge region 1 is the triangle
#' \eqn{T(B,CM,C)}, edge region 2 is \eqn{T(A,CM,C)} and edge region 3 is \eqn{T(A,B,CM)}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the edge regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with the elements
#' \item{re}{Indices (i.e., a \code{vector} of indices) of the edges whose region contains points in \code{Xp}
#' in the triangle \code{tri}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#' \item{desc}{Description of the edge labels as \code{"Edge labels are AB=3, BC=1, and AC=2"}.}
#'
#' @seealso \code{\link{rel.edges.tri.cent}}, \code{\link{rel.verts.tri.cent}} and \code{\link{rel.verts.tri.nd}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' P<-c(.4,.2)
#' rel.edges.triCM(P,Tr)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' re<-rel.edges.triCM(Xp,Tr)
#' re
#' CM<-(A+B+C)/3
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-Tr; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]+c(-.02,.03,.02)
#' yc<-Tr[,2]+c(.02,.02,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(CM,Ds)
#' xc<-txt[,1]+c(.05,.06,-.05,-.02)
#' yc<-txt[,2]+c(.03,.03,.05,-.08)
#' txt.str<-c("CM","re=2","re=3","re=1")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(re$re))
#' }
#'
#' @export rel.edges.triCM
rel.edges.triCM <- function(Xp,tri)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,]
  CM<-(A+B+C)/3;
  tri.ABM<-rbind(A,B,CM)
  tri.BCM<-rbind(B,C,CM)
  tri.ACM<-rbind(A,C,CM)

  nt<-nrow(Xp)
  if (nt>=1)
  {
    ind.set<-rep(NA,nt)
    for (i in 1:nt)
    {
      if (in.triangle(Xp[i,],tri.ABM,boundary = TRUE)$in.tri)
      {
        ind.set[i]<-3
      } else if (in.triangle(Xp[i,],tri.BCM,boundary = TRUE)$in.tri)
      {
        ind.set[i]<-1
      } else if (in.triangle(Xp[i,],tri.ACM,boundary = TRUE)$in.tri)
      {
        ind.set[i]<-2
      }
    }
  } else
  {
    {ind.set<-NA}
  }
  row.names(tri)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=ind.set, #relative edge
       tri=tri, #vertex labeling
       desc=edge.desc)
} #end of the function
#'

#################################################################

#' @title The indices of the \code{M}-edge regions in a triangle that contains the points in a give data set
#'
#' @description Returns the indices of the edges whose regions contain the points in data set \code{Xp} in
#' a triangle \code{tri}\eqn{=T(A,B,C)} and edge regions are based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}
#' (see the plots in the example for illustrations).
#'
#' The vertices of the  triangle \code{tri} are labeled as
#' \eqn{1=A}, \eqn{2=B}, and \eqn{3=C} also according to the row number the vertex is recorded in \code{tri}
#' and the corresponding edges are \eqn{1=BC}, \eqn{2=AC}, and \eqn{3=AB}.
#'
#' If a point in \code{Xp} is not inside \code{tri}, then the function yields \code{NA} as output for that entry.
#' The corresponding edge region is the polygon
#' with the vertex, \code{M}, and vertices other than the non-adjacent vertex, i.e., edge region 1 is the triangle
#' \eqn{T(B,M,C)}, edge region 2 is \eqn{T(A,M,C)} and edge region 3 is \eqn{T(A,B,M)}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the edge regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{re}{Indices (i.e., a \code{vector} of indices) of the edges whose region contains points in \code{Xp}
#' in the triangle \code{tri}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index opposite to edge
#' whose index is given in re.}
#' \item{desc}{Description of the edge labels as \code{"Edge labels are AB=3, BC=1, and AC=2"}.}
#'
#' @seealso \code{\link{rel.edges.triCM}}, \code{\link{rel.verts.tri.cent}} and \code{\link{rel.verts.tri.nd}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' M<-c(1.6,1.2)
#'
#' P<-c(.4,.2)
#' rel.edges.tri.cent(P,Tr,M)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' (re<-rel.edges.tri.cent(Xp,Tr,M))
#'
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Tr,pch=".",xlab="",ylab="",
#' main="Scatterplot of data points \n and the M-edge regions",axes=TRUE,
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-Tr; R<-rbind(M,M,M)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]+c(-.02,.03,.02)
#' yc<-Tr[,2]+c(.02,.02,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(.05,.06,-.05,-.02)
#' yc<-txt[,2]+c(.03,.03,.05,-.08)
#' txt.str<-c("M","re=2","re=3","re=1")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(re$re))
#' }
#'
#' @export rel.edges.tri.cent
rel.edges.tri.cent <- function(Xp,tri,M)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else  {
    Xp<-as.matrix(Xp)
    if (ncol(Xp)!=2 )
    {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,]
  tri.ABM<-rbind(A,B,M)
  tri.BCM<-rbind(B,C,M)
  tri.ACM<-rbind(A,C,M)

  nt<-nrow(Xp)
  if (nt>=1)
  {
    ind.set<-rep(NA,nt)
    for (i in 1:nt)
    {
      if (in.triangle(Xp[i,],tri.ABM,boundary = TRUE)$in.tri)
      {
        ind.set[i]<-3
      } else if (in.triangle(Xp[i,],tri.BCM,boundary = TRUE)$in.tri)
      {
        ind.set[i]<-1
      } else if (in.triangle(Xp[i,],tri.ACM,boundary = TRUE)$in.tri)
      {
        ind.set[i]<-2
      }
    }
  } else
  {
    {ind.set<-NA}
  }

  ind.set
  row.names(tri)<-c("A","B","C")  #vertex labeling
  edge.desc<-"Edge labels are AB=3, BC=1, and AC=2"

  list(re=ind.set, #relative edge
       tri=tri, #vertex labeling
       desc=edge.desc)
} #end of the function
#'

#################################################################

#' @title Region index inside the Gamma-1 region
#'
#' @description Returns the region index of the point \code{p} for the 6 regions in standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}, starting with 1 on the first one-sixth of the triangle,
#' and numbering follows the counter-clockwise direction (see the plot in the examples).
#' These regions are in the inner hexagon which is the Gamma-1 region for CS-PCD with \eqn{t=1}
#' if \code{p} is not in any of the 6 regions the function returns \code{NA}.
#'
#' @param p A 2D point whose index for the 6 regions in standard equilateral triangle \eqn{T_e} is determined.
#'
#' @return rel An integer between 1-6 (inclusive) or \code{NA}
#'
#' @seealso \code{\link{runif.std.tri.onesixth}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P<-c(.4,.2)
#' rel.six.Te(P)
#'
#' A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' CM<-(A+B+C)/3
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' h1<-c(1/2,sqrt(3)/18); h2<-c(2/3, sqrt(3)/9); h3<-c(2/3, 2*sqrt(3)/9);
#' h4<-c(1/2, 5*sqrt(3)/18); h5<-c(1/3, 2*sqrt(3)/9); h6<-c(1/3, sqrt(3)/9);
#'
#' r1<-(h1+h6+CM)/3;r2<-(h1+h2+CM)/3;r3<-(h2+h3+CM)/3;
#' r4<-(h3+h4+CM)/3;r5<-(h4+h5+CM)/3;r6<-(h5+h6+CM)/3;
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' L<-Te; R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' polygon(rbind(h1,h2,h3,h4,h5,h6))
#'
#' txt<-rbind(h1,h2,h3,h4,h5,h6)
#' xc<-txt[,1]+c(-.02,.02,.02,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,0,0,0)
#' txt.str<-c("h1","h2","h3","h4","h5","h6")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(Te,CM,r1,r2,r3,r4,r5,r6)
#' xc<-txt[,1]+c(-.02,.02,.02,0,0,0,0,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,0,0,0,0,0,0,0)
#' txt.str<-c("A","B","C","CM","1","2","3","4","5","6")
#' text(xc,yc,txt.str)
#'
#' n<-10  #try also n<-40
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' rsix<-vector()
#' for (i in 1:n)
#'   rsix<-c(rsix,rel.six.Te(Xp[i,]))
#' rsix
#'
#' plot(A,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp,pch=".")
#' L<-Te; R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' polygon(rbind(h1,h2,h3,h4,h5,h6))
#' text(Xp,labels=factor(rsix))
#'
#' txt<-rbind(Te,CM)
#' xc<-txt[,1]+c(-.02,.02,.02,0)
#' yc<-txt[,2]+c(.02,.02,.02,-.05)
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
rel.six.Te <- function(p)
{
  if (!is.point(p))
  {stop('The argument must be a numeric 2D point')}

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)
  if (in.triangle(p,Te,boundary = TRUE)$in.tri==FALSE)
  {rel<-NA; return(rel); stop}

  rel<-NA
  if (p[1] <= 1/2 && p[2] <= .5773502693*p[1] && p[2] >= 0.3849001795 - 0.5773502693*p[1])
    rel<-1
  else
  {
    if (p[1] > 1/2 && p[2] <= .5773502693-.5773502693*p[1] && p[2] >= -.1924500898+.5773502693*p[1])
      rel<-2
    else
    {
      if (p[1] <= 2/3 && p[2] > .5773502693-.5773502693*p[1] && p[2] <= .5773502693*p[1])
        rel<-3
      else
      {
        if (p[1] >= 1/2 && p[2] <= .7698003590-.5773502693*p[1] && p[2] > .5773502693*p[1])
          rel<-4
        else
        {
          if (p[1] < 1/2 && p[2] <= 0.1924500898 + 0.5773502693*p[1] && p[2] >= 0.5773502693 - 0.5773502693*p[1])
            rel<-5
          else
          {
            if (p[1] >= 1/3 && p[2] > .5773502693*p[1] && p[2] < 0.5773502693 - 0.5773502693*p[1])
              rel<-6
          }
        }
      }
    }
  }
  rel
} #end of the function
#'

#################################################################

# funsIndDelTri
#'
#' @title Functions provide the indices of the Delaunay triangles where the points reside
#'
#' @description
#' Two functions: \code{ind.del.tri} and \code{indices.del.tri}.
#'
#' \code{ind.del.tri} finds the index of the Delaunay triangle in which the given point, \code{p}, resides.
#'
#' \code{indices.del.tri} finds the indices of triangles for all the points in data set, \code{Xp}, as a vector.
#'
#' Delaunay triangulation is based on \code{Yp} and \code{DTmesh} are the Delaunay triangles with default \code{NULL}.
#' The function returns \code{NA} for a point not inside the convex hull of \code{Yp}.
#' Number of \code{Yp} points (i.e., size of \code{Yp}) should be at least three and the points
#' should be in general position so that Delaunay triangulation is (uniquely) defined.
#'
#' If the number of \code{Yp} points is 3, then there is only one Delaunay triangle and the indices of all
#' the points inside this triangle are all 1.
#'
#' See (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay
#' triangulation and the corresponding algorithm.
#'
#' @param p A 2D point; the index of the Delaunay triangle in which \code{p} resides is to be
#' determined. It is an argument for \code{ind.del.tri}.
#' @param Xp A set of 2D points representing the set of data points for which the indices of the Delaunay
#' triangles they reside is to be determined. It is an argument for \code{indices.del.tri}.
#' @param Yp A set of 2D points from which Delaunay triangulation is constructed.
#' @param DTmesh Delaunay triangles based on \code{Yp}, default is \code{NULL},
#' which is computed via \code{\link[interp]{tri.mesh}} function
#' in \code{interp} package. \code{\link[interp]{triangles}} function yields a triangulation data structure from the triangulation object
#' created by \code{tri.mesh}.
#'
#' @return \code{ind.del.tri} returns the index of the Delaunay triangle in which the given point, \code{p}, resides
#' and \code{indices.del.tri} returns the \code{vector} of indices of the Delaunay triangles in which points in the data
#' set, \code{Xp}, reside.
#'
#' @name funsIndDelTri
NULL
#'
#' @rdname funsIndDelTri
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for ind.del.tri
#' nx<-20 #number of X points (target)
#' ny<-5 #number of Y points (nontarget)
#' set.seed(1)
#' Yp<-cbind(runif(ny),runif(ny))
#'
#' Xp<-runifMT(nx,Yp)$g #data under CSR in the convex hull of Ypoints
#' #try also Xp<-cbind(runif(nx),runif(nx))
#'
#' ind.del.tri(Xp[10,],Yp)
#'
#' #or use
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation
#' TRY<-interp::triangles(DTY)[,1:3];
#' ind.del.tri(Xp[10,],Yp,DTY)
#'
#' ind.DT<-vector()
#' for (i in 1:nx)
#'  ind.DT<-c(ind.DT,ind.del.tri(Xp[i,],Yp))
#' ind.DT
#'
#' Xlim<-range(Yp[,1],Xp[,1])
#' Ylim<-range(Yp[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation based on Y points
#'
#' #plot of the data in the convex hull of Y points together with the Delaunay triangulation
#' plot(Xp,main=" ", xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),type="n")
#' interp::plot.triSht(DTY, add=TRUE, do.points = TRUE,pch=16,col="blue")
#' points(Xp,pch=".",cex=3)
#' text(Xp,labels = factor(ind.DT))
#' }
#'
#' @export
ind.del.tri <- function(p,Yp,DTmesh=NULL)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  ind<-NA;
  if (nrow(Yp)==3)
  {
    if (in.triangle(p,Yp)$in.tri)
    {ind<-1}
  } else
  {
    Yp<-as.matrix(Yp)
    if (!is.numeric(Yp) || ncol(Yp)!=2 || nrow(Yp)<3)
    {stop('Yp must be numeric and of dimension kx2 with k>=3')}

    if (is.null(DTmesh))
    {DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation
    }
    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3); #the Delaunay triangles

    nt<-nrow(DTr)  #number of Delaunay triangles;

    for (i in 1:nt)
    {
      tri<-Yp[DTr[i,],];
      if (in.triangle(p,tri)$in.tri)
      {ind<-i}
    }
  }
  ind
} #end of the function
#'
#' @rdname funsIndDelTri
#'
#' @examples
#' \dontrun{
#' #Examples for indices.del.tri
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Yp<-cbind(runif(ny),runif(ny))
#' Xp<-runifMT(nx,Yp)$g #data under CSR in the convex hull of Ypoints
#' #try also Xp<-cbind(runif(nx),runif(nx))
#'
#' tr.ind<-indices.del.tri(Xp,Yp)  #indices of the Delaunay triangles
#' tr.ind
#'
#' #or use
#' DTY<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation based on Y points
#' tr.ind<-indices.del.tri(Xp,Yp,DTY)  #indices of the Delaunay triangles
#' tr.ind
#'
#' Xlim<-range(Yp[,1],Xp[,1])
#' Ylim<-range(Yp[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' #plot of the data in the convex hull of Y points together with the Delaunay triangulation
#'
#' par(pty="s")
#' plot(Xp,main=" ", xlab=" ", ylab=" ",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),pch=".")
#' interp::plot.triSht(DTY, add=TRUE, do.points = TRUE,pch=16,col="blue")
#' text(Xp,labels = factor(tr.ind))
#' }
#'
#' @export
indices.del.tri <- function(Xp,Yp,DTmesh=NULL)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  nt<-nrow(Xp)

  Yp<-as.matrix(Yp)
  if (!is.numeric(Yp) || ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be numeric and of dimension kx2 with k>=3')}

 ind.set<-vector()
  for (i in 1:nt)
  {
   ind.set<-c(ind.set, ind.del.tri(Xp[i,],Yp,DTmesh) )
  }
 ind.set
} #end of the function
#'

#################################################################

#' @title The scatterplot of points from one class and plot of the Delaunay triangulation of the other class
#'
#' @description Plots the scatter plot of \code{Xp} points together with the Delaunay triangles based on the \code{Yp} points.
#' Both sets of points are of 2D.
#'
#' See (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points whose scatterplot is to be plotted.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges (default=\code{NULL} for both)
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A scatterplot of \code{Xp} points and the Delaunay triangulation of \code{Yp} points.
#'
#' @seealso \code{\link[interp]{plot.triSht}} in \code{interp} package
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-20; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' plotDeltri(Xp,Yp,xlab="",ylab="",main="X points and Delaunay Triangulation of Y points")
#' }
#'
#' @export plotDeltri
plotDeltri <- function(Xp,Yp,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('Xp and Yp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (is.null(xlim))
  {xlim<-range(Xp[,1],Yp[,1])}
  if (is.null(ylim))
  {ylim<-range(Xp[,2],Yp[,2])}

  oldpar <- par(no.readonly = TRUE)    # default par options
  on.exit(par(oldpar))            # sets default par options when the function is exited

  if (is.null(main))
  { main=paste(xname," points and Delaunay Triangulation of ",yname, " points",sep="")}

  if (nrow(Yp)==3)
  {
    par(mfrow=c(1,1),mar=c(5,5,4,2))
    plot(Xp[,1],Xp[,2],main=main, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=4,...)
    polygon(Yp,lty=2)
    } else
  {
    Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation
    par(mfrow=c(1,1),mar=c(5,5,4,2))
    plot(Xp[,1],Xp[,2],main=main, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=4,...)
    interp::plot.triSht(Ytrimesh, add=TRUE, do.points = TRUE,...)
    }
} #end of the function
#'

#################################################################

#' @title Points from one class inside the convex hull of the points from the other class
#'
#' @description Given two 2D data sets, \code{Xp} and \code{Yp}, it returns the \code{Xp} points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the data set.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#'
#' @return \code{Xp} points inside the convex hull of \code{Yp} points
#'
#' @seealso \code{\link{plotDeltri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' DT<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
#'
#' Xlim<-range(Xp[,1],Yp[,1])
#' Ylim<-range(Xp[,2],Yp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' Xch<-XinConvHullY(Xp,Yp)
#'
#' plot(Xp,main=" ", xlab=" ", ylab=" ",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),pch=".",cex=3)
#' interp::convex.hull(DT,plot.it = TRUE, add = TRUE)  # or try polygon(Yp[ch$i,])
#' points(Xch,pch=4,col="red")
#' }
#'
#' @export XinConvHullY
XinConvHullY <- function(Xp,Yp)
{
  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('both arguments must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  nx<-nrow(Xp)

  if (nrow(Yp)==3)
  {
    ch<-rep(0,nx)
    for (i in 1:nx)
    {ch[i]<-in.triangle(Xp[i,],Yp,boundary = TRUE)$in.tri}

    Xch<-Xp[ch==1,] #the Xp points inside the convex hull of Yp
  } else
  {
    DT<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    ch<-rep(0,nx)
    for (i in 1:nx)
    {ch[i]<-interp::in.convex.hull(DT,Xp[i,1],Xp[i,2],strict=FALSE)}

    Xch<-Xp[ch==1,] #the Xp points inside the convex hull of Yp
  }
  Xch
} #end of the function
#'

#################################################################

#' @title Number of Delaunay triangles based on a 2D data set
#'
#' @description Returns the number of Delaunay triangles based on the 2D set of points \code{Yp}.
#' See (\insertCite{okabe:2000,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and
#' the corresponding algorithm.
#'
#' @param Yp A set of 2D points which constitute the vertices of Delaunay triangles.
#'
#' @return Number of Delaunay triangles based on \code{Yp} points.
#'
#' @seealso \code{\link{plotDeltri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' ny<-10
#'
#' set.seed(1)
#' Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' num.del.tri(Yp)
#'
#' @export num.del.tri
num.del.tri <- function(Yp)
{
  if (!is.numeric(as.matrix(Yp)))
  {stop('the argument must be numeric')}

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('the argument must be of dimension kx2 with k>=3')}

  if (nrow(Yp)==3)
  {
    vec1<-rep(1,3);
    D0<-det(matrix(cbind(Yp,vec1),ncol=3))
    if (round(D0,14)==0)
    {stop('the Delaunay triangle is degenerate')}
    nt<-1
  } else
  {
    Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")  #Delaunay triangulation
    Ytri<-matrix(interp::triangles(Ytrimesh)[,1:3],ncol=3); #the Delaunay triangles
    nt<-nrow(Ytri)  #number of Delaunay triangles
  }
  nt
} #end of the function
#'

#################################################################

#' @title The indices of the \eqn{CM}-vertex regions in a triangle that contains the points in a give data set
#'
#' @description Returns the indices of the vertices whose regions contain the points in data set \code{Xp} in
#' a triangle \code{tri}\eqn{=(A,B,C)} and vertex regions are based on the center of mass \eqn{CM} of \code{tri}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the  triangle \code{tri} are labeled as
#' \eqn{1=A}, \eqn{2=B}, and \eqn{3=C} also according to the row number the vertex is recorded in \code{tri}.
#' If a point in \code{Xp} is not inside \code{tri}, then the function yields \code{NA} as output for that entry.
#' The corresponding vertex region is the polygon
#' with the vertex, \eqn{CM}, and midpoints the edges crossing the vertex.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the vertex regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Indices (i.e., a \code{vector} of indices) of the vertices whose region contains points in \code{Xp}
#' in the triangle \code{tri}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rel.verts.tri.cent}}, \code{\link{rel.verts.triCC}} and \code{\link{rel.verts.tri.nd}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' P<-c(.4,.2)
#' rel.verts.triCM(P,Tr)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' rv<-rel.verts.triCM(Xp,Tr)
#' rv
#'
#' CM<-(A+B+C)/3
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-matrix(rep(CM,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]+c(-.04,.05,.05)
#' yc<-Tr[,2]+c(-.05,.05,.03)
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(CM,Ds)
#' xc<-txt[,1]+c(.04,.04,-.03,0)
#' yc<-txt[,2]+c(-.07,.04,.06,-.08)
#' txt.str<-c("CM","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(rv$rv))
#' }
#'
#' @export rel.verts.triCM
rel.verts.triCM <- function(Xp,tri)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  nt<-nrow(Xp)
  if (nt==0)
  {ind.set<-NA } else
  {
    A<-tri[1,]; B<-tri[2,]; C<-tri[3,]
    CM<-(A+B+C)/3;
    M1<-(A+B)/2; M2<-(B+C)/2; M3<-(A+C)/2

    #t.mA<-interp::tri.mesh(c(A[1],M1[1],CM[1],M3[1]),c(A[2],M1[2],CM[2],M3[2]),duplicate="remove")
    #t.mB<-interp::tri.mesh(c(B[1],M2[1],CM[1],M1[1]),c(B[2],M2[2],CM[2],M1[2]),duplicate="remove")
    #t.mC<-interp::tri.mesh(c(C[1],M3[1],CM[1],M2[1]),c(C[2],M3[2],CM[2],M2[2]),duplicate="remove")

    #ind.vA<-interp::in.convex.hull(t.mA,Xp[,1],Xp[,2],strict=FALSE)
    #ind.vB<-interp::in.convex.hull(t.mB,Xp[,1],Xp[,2],strict=FALSE)
    #ind.vC<-interp::in.convex.hull(t.mC,Xp[,1],Xp[,2],strict=FALSE)

    if (length(unique(Xp))==length(Xp))
    {
    ind.vA<-gMOIP::inHull(Xp, rbind(A,M1,CM,M3))>=0
    ind.vB<-gMOIP::inHull(Xp, rbind(B,M2,CM,M1))>=0
    ind.vC<-gMOIP::inHull(Xp, rbind(C,M3,CM,M2))>=0
    } else
    {
    ind.vA<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(A,M1,CM,M3)))>=0
    ind.vB<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(B,M2,CM,M1)))>=0
    ind.vC<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(C,M3,CM,M2)))>=0
    }

    ind.set<-rep(NA,nt)
    ind.set[ind.vA==TRUE]<-1
    ind.set[ind.vB==TRUE]<-2
    ind.set[ind.vC==TRUE]<-3
  }

  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=ind.set, #relative vertices
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the vertex region in a standard basic triangle form that contains a given point
#'
#' @description Returns the index of the related vertex in the standard basic triangle form whose region contains point \code{p}.
#' The standard basic triangle form is \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}..
#'
#' Vertex regions are based on the general center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the standard basic triangle form \eqn{T_b}.
#' Vertices of the standard basic triangle form \eqn{T_b} are labeled according to the row number the
#' vertex is recorded, i.e., as 1=(0,0), 2=(1,0),and \eqn{3=(c_1,c_2)}.
#'
#' If the point, \code{p}, is not inside \eqn{T_b}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon with the vertex, \code{M}, and projections from \code{M}
#' to the edges on the lines joining vertices and \code{M}.
#' That is, \code{rv=1} has vertices \eqn{(0,0),D_3,M,D_2}; \code{rv=2} has vertices \eqn{(1,0),D_1,M,D_3}; and
#' \eqn{rv=3} has vertices \eqn{(c_1,c_2),D_2,M,D_1} (see the illustration in the examples).
#'
#' Any given triangle can be mapped to the standard basic triangle form
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \code{M}-vertex region it resides in is to be determined in the
#' standard basic triangle form \eqn{T_b}.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle form
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle form.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the vertex whose region contains point, \code{p}; index of the vertex is the same as the row
#' number in the standard basic triangle form, \eqn{T_b}}
#' \item{tri}{The vertices of the standard basic triangle form, \eqn{T_b}, where row number corresponds to the vertex index \code{rv}
#' with \code{rv=1} is row \eqn{1=(0,0)}, \code{rv=2} is row \eqn{2=(1,0)}, and \eqn{rv=3} is row \eqn{3=(c_1,c_2)}.}
#'
#' @seealso \code{\link{rv.bas.triCM}}, \code{\link{rv.tri.cent}}, \code{\link{rv.triCC}},
#' \code{\link{rv.bas.triCC}}, \code{\link{rv.triCM}}, and \code{\link{rvTeCM}}
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
#' M<-c(.6,.2)
#'
#' P<-c(.4,.2)
#' rv.bas.tri.cent(P,c1,c2,M)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' Rv<-vector()
#' for (i in 1:n)
#' { Rv<-c(Rv,rv.bas.tri.cent(Xp[i,],c1,c2,M)$rv)}
#' Rv
#'
#' Ds<-cent2edges.bas.tri(c1,c2,M)
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tb)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Tb,pch=".",xlab="",ylab="",axes=TRUE,
#' xlim=Xlim+xd*c(-.1,.1),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' points(Xp,pch=".",col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tb[,1]+c(-.04,.05,.04)
#' yc<-Tb[,2]+c(.02,.02,.03)
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(-.02,.04,-.03,0)
#' yc<-txt[,2]+c(-.02,.02,.02,-.03)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(Rv))
#' }
#'
#' @export rv.bas.tri.cent
rv.bas.tri.cent <- function(p,c1,c2,M)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (in.triangle(M,Tb,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the standard basic triangle form')}

  a1<-y1[1]; a2<-y1[2]; b1<-y2[1]; b2<-y2[2]; c1<-y3[1]; c2<-y3[2];

  if (in.triangle(p,Tb)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    Ms<-cent2edges.bas.tri(c1,c2,M)
    M1<-Ms[1,]; M2<-Ms[2,]; M3<-Ms[3,];

    if (in.triangle(p,rbind(y1,M3,M),boundary=TRUE)$in.tri | in.triangle(p,rbind(y1,M,M2),boundary=TRUE)$in.tri)
    {rv<-1}
    else
    {
      if (in.triangle(p,rbind(M3,y2,M),boundary=TRUE)$in.tri | in.triangle(p,rbind(y2,M1,M),boundary=TRUE)$in.tri)
      {rv<-2}
      else
      {rv<-3}
    }
  }
  row.names(Tb)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=Tb #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CC}-vertex region in a standard basic triangle form that contains a point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' the standard basic triangle form \eqn{T_b=T((0,0),(1,0),(c_1,c_2))} where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}
#' and vertex regions are based on the circumcenter \eqn{CC} of \eqn{T_b}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the standard basic triangle form \eqn{T_b} are labeled as
#' \eqn{1=(0,0)}, \eqn{2=(1,0)},and \eqn{3=(c_1,c_2)} also according to the row number the vertex is recorded in \eqn{T_b}.
#' If the point, \code{p}, is not inside \eqn{T_b}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon whose interior points are closest to that vertex.
#'
#' Any given triangle can be mapped to the standard basic triangle form
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle form is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CC}-vertex region it resides in is to be determined in the
#' standard basic triangle form \eqn{T_b}.
#' @param c1,c2 Positive real numbers which constitute the upper vertex of the standard basic triangle form (i.e., the vertex
#' adjacent to the shorter edges of \eqn{T_b}); \eqn{c_1} must be in \eqn{[0,1/2]},
#' \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the \eqn{CC}-vertex region that contains point, \code{p} in the standard basic triangle form \eqn{T_b}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}
#' with row \eqn{1=(0,0)}, row \eqn{2=(1,0)}, and row \eqn{3=(c_1,c_2)}.}
#'
#' @seealso \code{\link{rv.triCM}}, \code{\link{rv.tri.cent}}, \code{\link{rv.triCC}},
#' \code{\link{rv.bas.triCM}}, \code{\link{rv.bas.tri.cent}}, and \code{\link{rvTeCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6;  #try also c1<-.5; c2<-.5;
#'
#' P<-c(.3,.2)
#' rv.bas.triCC(P,c1,c2)
#'
#' A<-c(0,0);B<-c(1,0);C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tb[,1])
#' Ylim<-range(Tb[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,asp=1,xlab="",ylab="",axes=TRUE,pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tb,CC,Ds)
#' xc<-txt[,1]+c(-.03,.03,0.02,-.01,.05,-.05,.01)
#' yc<-txt[,2]+c(.02,.02,.03,.06,.03,.03,-.03)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' RV1<-(A+D3+CC+D2)/4
#' RV2<-(B+D3+CC+D1)/4
#' RV3<-(C+D2+CC+D1)/4
#'
#' txt<-rbind(RV1,RV2,RV3)
#' xc<-txt[,1]
#' yc<-txt[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rv.bas.triCC(Xp[i,],c1,c2)$rv)
#' Rv
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,asp=1,xlab="",pch=".",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' points(Xp,pch=".")
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(Rv))
#'
#' txt<-rbind(Tb,CC,Ds)
#' xc<-txt[,1]+c(-.03,.03,0.02,-.01,.05,-.05,.01)
#' yc<-txt[,2]+c(.02,.02,.03,.06,.03,.03,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export rv.bas.triCC
rv.bas.triCC <- function(p,c1,c2)
{
  if (!is.point(p))
  {stop('p must be a numeric point of dimension 2')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  p1<-c(0,0); p2<-c(1,0); p3<-c(c1,c2);
  Tb<-rbind(p1,p2,p3)
  if (in.triangle(p,Tb,boundary=TRUE)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    x<-p[1]; y<-p[2];
    if ( y<= lD2CCinTb(x,c1,c2)$y && x<=1/2)
    {rv<-1}
    else
    {
      if (y<= lD1CCinTb(x,c1,c2)$y && x>=1/2)
      {rv<-2}
      else
      {rv<-3}
    }
  }
  row.names(Tb)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=Tb #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the vertex region in a triangle that contains a given point
#'
#' @description Returns the index of the related vertex in the triangle, \code{tri}, whose region contains point \code{p}.
#'
#' Vertex regions are based on the general center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}.
#' Vertices of the triangle \code{tri} are labeled according to the row number the vertex is recorded.
#'
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon with the vertex, \code{M}, and projections from \code{M}
#' to the edges on the lines joining vertices and \code{M}. (see the illustration in the examples).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \code{M}-vertex region it resides in is to be determined in the
#' triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the vertex whose region contains point, \code{p}; index of the vertex is the same as the row
#' number in the triangle, \code{tri}}
#' \item{tri}{The vertices of the triangle, \code{tri}, where row number corresponds to the vertex index \code{rv}
#' with \code{rv=1} is row 1, \code{rv=2} is row 2, and \eqn{rv=3} is is row 3.}
#'
#' @seealso \code{\link{rv.triCM}}, \code{\link{rv.triCC}}, \code{\link{rv.bas.triCC}},
#' \code{\link{rv.bas.triCM}}, \code{\link{rv.bas.tri.cent}}, and \code{\link{rvTeCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' M<-c(1.6,1.0)
#'
#' P<-c(1.5,1.6)
#' rv.tri.cent(P,Tr,M)
#' #try also rv.tri.cent(P,Tr,M=c(2,2))  #center is not in the interior of the triangle
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' Rv<-vector()
#' for (i in 1:n)
#' {Rv<-c(Rv,rv.tri.cent(Xp[i,],Tr,M)$rv)}
#' Rv
#'
#' Ds<-cent2edges.tri(Tr,M)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Tr,pch=".",xlab="",ylab="",main="Illustration of M-Vertex Regions\n in a Triangle",axes=TRUE,
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]
#' yc<-Tr[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(-.02,.04,-.04,0)
#' yc<-txt[,2]+c(-.02,.04,.05,-.08)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(Rv))
#' }
#'
#' @export rv.tri.cent
rv.tri.cent <- function(p,tri,M)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];
  a1<-y1[1]; a2<-y1[2]; b1<-y2[1]; b2<-y2[2]; c1<-y3[1]; c2<-y3[2];

  if (in.triangle(p,tri,boundary = T)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
    {stop('center is not in the interior of the triangle')}

    Ds<-cent2edges.tri(tri,M)
    D1<-Ds[1,]; D2<-Ds[2,]; D3<-Ds[3,];

    if (in.triangle(p,rbind(y1,D3,M),boundary=TRUE)$in.tri | in.triangle(p,rbind(y1,M,D2),boundary=TRUE)$in.tri)
    {rv<-1}
    else
    {
      if (in.triangle(p,rbind(D3,y2,M),boundary=TRUE)$in.tri | in.triangle(p,rbind(y2,D1,M),boundary=TRUE)$in.tri)
      {rv<-2}
      else
      {rv<-3}
    }
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CC}-vertex region in a triangle that contains a point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' a triangle \code{tri}\eqn{=(A,B,C)} and vertex regions are based on the circumcenter \eqn{CC} of \code{tri}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the  triangle \code{tri} are labeled as
#' \eqn{1=A}, \eqn{2=B}, and \eqn{3=C} also according to the row number the vertex is recorded in \code{tri}.
#' If the point, \code{p}, is not inside \code{tri}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon whose interior points are closest to that vertex.
#' If \code{tri} is equilateral triangle, then \eqn{CC} and \eqn{CM} (center of mass) coincide.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CC}-vertex region it resides in is to be determined in the
#' triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the \eqn{CC}-vertex region that contains point, \code{p} in the triangle \code{tri}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rv.tri.cent}}, \code{\link{rv.triCM}}, \code{\link{rv.bas.triCM}},
#' \code{\link{rv.bas.triCC}}, \code{\link{rv.bas.tri.cent}}, and \code{\link{rvTeCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' P<-c(1.3,1.2)
#' rv.triCC(P,Tr)
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],CC[1])
#' Ylim<-range(Tr[,2],CC[2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,asp=1,xlab="",ylab="",pch=".",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tr,CC,Ds)
#' xc<-txt[,1]+c(-.07,.08,.06,.12,-.1,-.1,-.09)
#' yc<-txt[,2]+c(.02,-.02,.03,.0,.02,.06,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' RV1<-(A+.5*(D3-A)+A+.5*(D2-A))/2
#' RV2<-(B+.5*(D3-B)+B+.5*(D1-B))/2
#' RV3<-(C+.5*(D2-C)+C+.5*(D1-C))/2
#'
#' txt<-rbind(RV1,RV2,RV3)
#' xc<-txt[,1]
#' yc<-txt[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' n<-20  #try also n<-40
#' Xp<-runif.tri(n,Tr)$g
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rv.triCC(Xp[i,],Tr)$rv)
#' Rv
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,asp=1,xlab="",ylab="",main="Illustration of CC-Vertex Regions\n in a Triangle",
#' pch=".",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".")
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' text(Xp,labels=factor(Rv))
#'
#' txt<-rbind(Tr,CC,Ds)
#' xc<-txt[,1]+c(-.07,.08,.06,.12,-.1,-.1,-.09)
#' yc<-txt[,2]+c(.02,-.02,.03,.0,.02,.06,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export rv.triCC
rv.triCC <- function(p,tri)
{
  if (!is.point(p))
  {stop('p must be a numeric point of dimension 2')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  p1<-tri[1,]; p2<-tri[2,]; p3<-tri[3,]
  if (in.triangle(p,tri,boundary=TRUE)$in.tri==FALSE)
  {rv<-NA
  } else
  {  mdt<-max(Dist(p1,p2),Dist(p1,p3),Dist(p2,p3)); #max edge length
  for (i in 1:3)
  {
    d1<-Dist(p,tri[i,]);
    if (d1<mdt)
    {mdt<-d1; rv<-i }
  }
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the vertex region in the standard equilateral triangle that contains a given point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with vertex regions are constructed with center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of triangle, \eqn{T_e}, are labeled as \eqn{1,2,3}
#' according to the row number the vertex is recorded in \eqn{T_e}. If the point, \code{p}, is not inside \eqn{T_e}, then the
#' function yields \code{NA} as output. The corresponding vertex region is the polygon with the vertex, \code{M}, and
#' projections from \code{M} to the edges on the lines joining vertices and \code{M}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \code{M}-vertex region it resides in is to be determined in the
#' standard equilateral triangle \eqn{T_e}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the vertex whose region contains point, \code{p}.}
#' \item{tri}{The vertices of the triangle, \eqn{T_e}, where row number corresponds to the vertex index in \code{rv}
#' with row \eqn{1=(0,0)}, row \eqn{2=(1,0)}, and row \eqn{3=(1/2,\sqrt{3}/2)}.}
#'
#' @seealso \code{\link{rvTeCM}}, \code{\link{rv.tri.cent}}, \code{\link{rv.triCC}},
#' \code{\link{rv.bas.triCC}}, \code{\link{rv.triCM}}, and \code{\link{rv.bas.tri.cent}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' n<-20  #try also n<-40
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' rvTe.cent(Xp[1,],M)
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rvTe.cent(Xp[i,],M)$rv)
#' Rv
#'
#' Ds<-cent2edges.tri(Te,M)
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Te)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",axes=TRUE,
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp,pch=".",col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Te,M)
#' xc<-txt[,1]+c(-.02,.03,.02,0)
#' yc<-txt[,2]+c(.02,.02,.03,.05)
#' txt.str<-c("A","B","C","M")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(Rv))
#' }
#'
#' @export rvTe.cent
rvTe.cent <- function(p,M)
{
  if (!is.point(p))
  {stop('the argument must be a numeric 2D point')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(1/2,sqrt(3)/2); Te<-rbind(y1,y2,y3);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(p,Te)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
    {stop('center is not in the interior of the triangle')}

    Ds<-cent2edges.tri(Te,M)
    D1<-Ds[1,]; D2<-Ds[2,]; D3<-Ds[3,];

    if (in.triangle(p,rbind(y1,D3,M),boundary=TRUE)$in.tri | in.triangle(p,rbind(y1,M,D2),boundary=TRUE)$in.tri)
    {rv<-1}
    else
    {
      if (in.triangle(p,rbind(D3,y2,M),boundary=TRUE)$in.tri | in.triangle(p,rbind(y2,D1,M),boundary=TRUE)$in.tri)
      {rv<-2}
      else
      {rv<-3}
    }
  }
  row.names(Te)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=Te #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CM}-vertex region in the standard equilateral triangle that contains a given point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with vertex regions are constructed with center of mass CM
#' (see the plots in the example for illustrations).
#'
#' The vertices of triangle, \eqn{T_e}, are labeled as \eqn{1,2,3}
#' according to the row number the vertex is recorded in \eqn{T_e}. If the point, \code{p}, is not inside \eqn{T_e}, then the
#' function yields \code{NA} as output. The corresponding vertex region is the polygon with the vertex, \eqn{CM}, and
#' midpoints of the edges adjacent to the vertex.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point for which \eqn{CM}-vertex region it resides in is to be determined in the
#' standard equilateral triangle \eqn{T_e}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the vertex whose region contains point, \code{p}.}
#' \item{tri}{The vertices of the triangle, \eqn{T_e}, where row number corresponds to the vertex index in \code{rv}
#' with row \eqn{1=(0,0)}, row \eqn{2=(1,0)}, and row \eqn{3=(1/2,\sqrt{3}/2)}.}
#'
#' @seealso \code{\link{rv.bas.triCM}}, \code{\link{rv.tri.cent}}, \code{\link{rv.triCC}},
#' \code{\link{rv.bas.triCC}}, \code{\link{rv.triCM}}, and \code{\link{rv.bas.tri.cent}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#'
#' n<-20  #try also n<-40
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' rvTeCM(Xp[1,])
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rvTeCM(Xp[i,])$rv)
#' Rv
#'
#' CM<-(A+B+C)/3
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,asp=1,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp,pch=".",col=1)
#' L<-matrix(rep(CM,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Te,CM)
#' xc<-txt[,1]+c(-.02,.03,.02,0)
#' yc<-txt[,2]+c(.02,.02,.03,.05)
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(Rv))
#' }
#'
#' @export rvTeCM
rvTeCM <- function(p)
{
  if (!is.point(p))
  {stop('the argument must be a numeric 2D point')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(1/2,sqrt(3)/2);
  tri<-rbind(y1,y2,y3);
  if (in.triangle(p,tri,boundary = TRUE)$in.tri==FALSE)
  {rv<-NA
  } else
  {
    mdt<-1; #maximum distance from a point inside Te to a vertex
    for (i in 1:3)
    {
      d1<-Dist(p,tri[i,]);
      if (d1<mdt)
      {mdt<-d1; rv<-i }
    }
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertex
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The indices of the vertex regions in a triangle that contains the points in a give data set
#'
#' @description Returns the indices of the vertices whose regions contain the points in data set \code{Xp} in
#' a triangle \code{tri}\eqn{=T(A,B,C)}.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)} in Cartesian
#' coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle
#' to the edges on the extension of the lines joining \code{M} to the vertices or based on the circumcenter of \code{tri}.
#' Vertices of triangle \code{tri} are labeled as \eqn{1,2,3} according to the row number the vertex is recorded.
#'
#' If a point in \code{Xp} is not inside \code{tri}, then the function yields \code{NA} as output for that entry.
#' The corresponding vertex region is the polygon with the vertex, \code{M}, and
#' projection points from \code{M} to the edges crossing the vertex (as the output of \code{cent2edges.tri(Tr,M)})
#' or \eqn{CC}-vertex region. (see the examples for an illustration).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:dom-num-NPE-Spat2011,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the vertex regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri} or the circumcenter of \code{tri}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Indices of the vertices whose regions contains points in \code{Xp}.}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rel.verts.triCM}}, \code{\link{rel.verts.triCC}} and \code{\link{rel.verts.tri.nd}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' M<-c(1.6,1.0)
#'
#' P<-c(.4,.2)
#' rel.verts.tri.cent(P,Tr,M)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also  #M<-c(1.6,1.0)
#'
#' rel.verts.tri.cent(Xp,Tr,M)
#' rel.verts.tri.cent(rbind(Xp,c(2,2)),Tr,M)
#'
#' rv<-rel.verts.tri.cent(Xp,Tr,M)
#' rv
#'
#' ifelse(identical(M,circ.cent.tri(Tr)),
#' Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2),Ds<-cent2edges.tri(Tr,M))
#'
#' Xlim<-range(Tr[,1],M[1],Xp[,1])
#' Ylim<-range(Tr[,2],M[2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' plot(Tr,pch=".",xlab="",ylab="",
#' main="Scatterplot of data points \n and M-vertex regions in a triangle",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]
#' yc<-Tr[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(.02,.04,-.03,0)
#' yc<-txt[,2]+c(.07,.04,.05,-.07)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(rv$rv))
#' }
#'
#' @export rel.verts.tri.cent
rel.verts.tri.cent <- function(Xp,tri,M)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  CC = circ.cent.tri(tri)
  if ( !(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri) )
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  nt<-nrow(Xp)

  if (nt==0)
  {ind.set<-NA}
  else
  {
    if (isTRUE(all.equal(M,CC)))
    {
      ind.set<-rel.verts.triCC(Xp,tri)$rv
    } else
    {
      A<-tri[1,]; a1<-A[1]; a2<-A[2];
      B<-tri[2,]; b1<-B[1]; b2<-B[2];
      C<-tri[3,]; c1<-C[1]; c2<-C[2];

      m1<-M[1]; m2<-M[2]

      M11<--(a1*b1*c2-a1*b1*m2-a1*b2*c1+a1*c1*m2+a2*b1*m1-a2*c1*m1-b1*c2*m1+b2*c1*m1)/(b2*a1-c2*a1-a2*b1+a2*c1+m2*b1-b2*m1-m2*c1+c2*m1);
      M12<-(a1*b2*m2-a1*c2*m2-a2*b1*c2+a2*b2*c1-a2*b2*m1+a2*c2*m1+b1*c2*m2-b2*c1*m2)/(b2*a1-c2*a1-a2*b1+a2*c1+m2*b1-b2*m1-m2*c1+c2*m1);
      D1<-c(M11,M12)
      M21 <-(a1*b1*c2-a1*b1*m2+a1*b2*m1-a1*c2*m1-a2*b1*c1+a2*c1*m1+b1*c1*m2-b2*c1*m1)/(b2*a1-m2*a1-a2*b1+m1*a2+c2*b1-c1*b2+m2*c1-c2*m1);
      M22 <-(a1*b2*c2-a1*c2*m2-a2*b1*m2-a2*b2*c1+a2*b2*m1+a2*c1*m2+b1*c2*m2-b2*c2*m1)/(b2*a1-m2*a1-a2*b1+m1*a2+c2*b1-c1*b2+m2*c1-c2*m1);
      D2<-c(M21,M22)
      M31 <- (a1*b2*c1-a1*b2*m1-a1*c1*m2+a1*c2*m1-a2*b1*c1+a2*b1*m1+b1*c1*m2-b1*c2*m1)/(c2*a1-m2*a1-a2*c1+m1*a2-c2*b1+m2*b1+c1*b2-b2*m1);
      M32 <- (a1*b2*c2-a1*b2*m2-a2*b1*c2+a2*b1*m2-a2*c1*m2+a2*c2*m1+b2*c1*m2-b2*c2*m1)/(c2*a1-m2*a1-a2*c1+m1*a2-c2*b1+m2*b1+c1*b2-b2*m1);
      D3<-c(M31,M32)

      # t.mA<-interp::tri.mesh(c(A[1],D3[1],M[1],D2[1]),c(A[2],D3[2],M[2],D2[2]),duplicate="error")
      # t.mB<-interp::tri.mesh(c(B[1],D1[1],M[1],D3[1]),c(B[2],D1[2],M[2],D3[2]),duplicate="error")
      # t.mC<-interp::tri.mesh(c(C[1],D2[1],M[1],D1[1]),c(C[2],D2[2],M[2],D1[2]),duplicate="error")
      #
      # ind.set<-rep(NA,nt)
      # ind.vA<-interp::in.convex.hull(t.mA,Xp[,1],Xp[,2], strict=FALSE)
      # ind.vB<-interp::in.convex.hull(t.mB,Xp[,1],Xp[,2], strict=FALSE)
      # ind.vC<-interp::in.convex.hull(t.mC,Xp[,1],Xp[,2], strict=FALSE)

      if (length(unique(Xp))==length(Xp))
      {
        ind.vA<-gMOIP::inHull(Xp, rbind(A,D3,M,D2))>=0
        ind.vB<-gMOIP::inHull(Xp, rbind(B,D1,M,D3))>=0
        ind.vC<-gMOIP::inHull(Xp, rbind(C,D2,M,D1))>=0
      } else
      {
        ind.vA<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(A,D3,M,D2)))>=0
        ind.vB<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(B,D1,M,D3)))>=0
        ind.vC<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(C,D2,M,D1)))>=0
      }

      ind.set<-rep(NA,nt)
      ind.set[ind.vA==TRUE]<-1
      ind.set[ind.vB==TRUE]<-2
      ind.set[ind.vC==TRUE]<-3
    }
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=ind.set, #relative vertices
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The alternative function for the indices of the M-vertex regions in a triangle that contains the points
#' in a give data set
#'
#' @description An alternative function to the function \code{\link{rel.verts.tri.cent}}
#' when the center M is not the circumcenter falling outside the triangle.
#' This function only works for a center \eqn{M} in the interior of the triangle,
#' with the projections of \eqn{M} to the edges along the lines joining \eqn{M} to the vertices.
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the vertex regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Indices of the vertices whose regions contains points in \code{Xp}.}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rel.verts.tri.cent}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' M<-c(1.6,1.0)
#'
#' P<-c(.4,.2)
#' rel.verts.triM(P,Tr,M)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-c(1.6,1.0)  #try also M<-c(1.3,1.3)
#'
#' (rv<-rel.verts.triM(Xp,Tr,M))
#' rel.verts.triM(rbind(Xp,c(2,2)),Tr,M)
#'
#' Ds<-cent2edges.tri(Tr,M)
#'
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]+c(-.03,.05,.05)
#' yc<-Tr[,2]+c(-.06,.02,.05)
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(.02,.04,-.03,0)
#' yc<-txt[,2]+c(.07,.03,.05,-.07)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(rv$rv))
#' }
#'
#' @export
rel.verts.triM <- function(Xp,tri,M)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the standard basic triangle form')}

  nt<-nrow(Xp)

  if (nt==0)
  {ind.set<-NA}
  else
  {
    rv<-rep(NA,nt)
    for (i in 1:nt)
    { if (in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri)
      rv[i]<-rv.tri.cent(Xp[i,],tri,M)$rv
    }
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=rv, #relative vertices
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The indices of the vertex regions in a triangle that contains the points in a give data set
#'
#' @description Returns the indices of the vertices whose regions contain the points in data set \code{Xp} in
#' a triangle \code{tri}\eqn{=(A,B,C)} and vertex regions are based on the center \code{cent} which yields nondegenerate asymptotic
#' distribution of the domination number of PE-PCD for uniform data in \code{tri} for expansion parameter \code{r} in \eqn{(1,1.5]}.
#'
#' Vertices of triangle \code{tri} are labeled as \eqn{1,2,3} according to the row number the vertex is recorded
#' if a point in \code{Xp} is not inside \code{tri}, then the function yields \code{NA} as output for that entry.
#' The corresponding vertex region is the polygon with the vertex, \code{cent}, and projection points on the edges.
#' The center label \code{cent} values \code{1,2,3} correspond to the vertices \eqn{M_1}, \eqn{M_2}, and \eqn{M_3}; with default 1
#' (see the examples for an illustration).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:dom-num-NPE-Spat2011,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the vertex regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]} for this function.
#' @param cent Index of the center (as \eqn{1,2,3} corresponding to \eqn{M_1,\,M_2,\,M_3}) which gives nondegenerate asymptotic
#' distribution of the domination number of PE-PCD for uniform data in \code{tri} for expansion parameter \code{r} in \eqn{(1,1.5]};
#' default \code{cent=1}.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Indices (i.e., a \code{vector} of indices) of the vertices whose region contains points in \code{Xp}
#' in the triangle \code{tri}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rel.verts.triCM}}, \code{\link{rel.verts.triCC}} and \code{\link{rel.verts.tri.cent}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' r<-1.35
#' cent<-2
#'
#' P<-c(1.4,1.0)
#' rel.verts.tri.nd(P,Tr,r,cent)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' rel.verts.tri.nd(Xp,Tr,r,cent)
#' rel.verts.tri.nd(rbind(Xp,c(2,2)),Tr,r,cent)
#'
#' rv<-rel.verts.tri.nd(Xp,Tr,r,cent)
#'
#' M<-cent.nondegPE(Tr,r)[cent,];
#' Ds<-cent.ndPE2edges(Tr,r,cent)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]+c(-.03,.05,.05)
#' yc<-Tr[,2]+c(-.06,.02,.05)
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(.02,.04,-.03,0)
#' yc<-txt[,2]+c(.07,.03,.05,-.07)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(rv$rv))
#' }
#'
#' @export rel.verts.tri.nd
rel.verts.tri.nd <- function(Xp,tri,r,cent=1)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(r,1) || r<=1 || r>1.5)
  {stop('r must be a scalar in (1,1.5]')}

  if (cent!=1 & cent!=2 & cent!=3)
  {stop('center index, cent, must be 1, 2 or 3')}

  Xp<-matrix(Xp,ncol=2)
  nt<-nrow(Xp)

  if (nt==0)
  {ind.set<-NA}
  else
  {
    A<-tri[1,]; B<-tri[2,]; C<-tri[3,];
    M<-cent.nondegPE(tri,r)[cent,];
    Ds<-cent.ndPE2edges(tri,r,cent)
    D1<-Ds[1,]; D2<-Ds[2,]; D3<-Ds[3,];

    # t.mA<-interp::tri.mesh(c(A[1],D3[1],M[1],D2[1]),c(A[2],D3[2],M[2],D2[2]),duplicate="remove")
    # t.mB<-interp::tri.mesh(c(B[1],D1[1],M[1],D3[1]),c(B[2],D1[2],M[2],D3[2]),duplicate="remove")
    # t.mC<-interp::tri.mesh(c(C[1],D2[1],M[1],D1[1]),c(C[2],D2[2],M[2],D1[2]),duplicate="remove")
    # ind.set<-rep(NA,nt)
    # ind.vA<-interp::in.convex.hull(t.mA,Xp[,1],Xp[,2],strict=FALSE)
    # ind.vB<-interp::in.convex.hull(t.mB,Xp[,1],Xp[,2],strict=FALSE)
    # ind.vC<-interp::in.convex.hull(t.mC,Xp[,1],Xp[,2],strict=FALSE)

    if (length(unique(Xp))==length(Xp))
    {
      ind.vA<-gMOIP::inHull(Xp, rbind(A,D3,M,D2))>=0
      ind.vB<-gMOIP::inHull(Xp, rbind(B,D1,M,D3))>=0
      ind.vC<-gMOIP::inHull(Xp, rbind(C,D2,M,D1))>=0
    } else
    {
      ind.vA<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(A,D3,M,D2)))>=0
      ind.vB<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(B,D1,M,D3)))>=0
      ind.vC<-apply(Xp,1, function(x) gMOIP::inHull(x,vertices=rbind(C,D2,M,D1)))>=0
    }

    ind.set<-rep(NA,nt)
    ind.set[ind.vA==TRUE]<-1
    ind.set[ind.vB==TRUE]<-2
    ind.set[ind.vC==TRUE]<-3
  }
  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=ind.set, #relative vertices
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The indices of the \eqn{CC}-vertex regions in a triangle that contains the points in a give data set.
#'
#' @description Returns the indices of the vertices whose regions contain the points in data set \code{Xp} in
#' a triangle \code{tri}\eqn{=(A,B,C)} and vertex regions are based on the circumcenter \eqn{CC} of \code{tri}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the  triangle \code{tri} are labeled as
#' \eqn{1=A}, \eqn{2=B}, and \eqn{3=C} also according to the row number the vertex is recorded in \code{tri}.
#' If a point in \code{Xp} is not inside \code{tri}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon whose interior points are closest to that vertex.
#' If \code{tri} is equilateral triangle, then \eqn{CC} and \eqn{CM} (center of mass) coincide.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points for which indices of the vertex regions
#' containing them are to be determined.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Indices (i.e., a \code{vector} of indices) of the vertices whose region contains points in \code{Xp}
#' in the triangle \code{tri}}
#' \item{tri}{The vertices of the triangle, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rel.verts.triCM}}, \code{\link{rel.verts.tri.cent}}, and \code{\link{rel.verts.tri.nd}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' P<-c(.4,.2)
#' rel.verts.triCC(P,Tr)
#'
#' n<-20  #try also n<-40
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' rel.verts.triCC(Xp,Tr)
#' rel.verts.triCC(rbind(Xp,c(2,2)),Tr)
#'
#' (rv<-rel.verts.triCC(Xp,Tr))
#'
#' CC<-circ.cent.tri(Tr)
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1],CC[1])
#' Ylim<-range(Tr[,2],Xp[,2],CC[2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",asp=1,xlab="",ylab="",
#' main="Scatterplot of data points \n and the CC-vertex regions",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=".",col=1)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' xc<-Tr[,1]
#' yc<-Tr[,2]
#' txt.str<-c("rv=1","rv=2","rv=3")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(CC,Ds)
#' xc<-txt[,1]+c(.04,.04,-.03,0)
#' yc<-txt[,2]+c(-.07,.04,.06,-.08)
#' txt.str<-c("CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' text(Xp,labels=factor(rv$rv))
#' }
#'
#' @export rel.verts.triCC
rel.verts.triCC <- function(Xp,tri)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  nt<-nrow(Xp)
  ind.set<-c()
  if (nt==0)
  {ind.set<-NA}
  else
  {
    for (i in 1:nt)
    {
      ind.set<-c(ind.set,rv.triCC(Xp[i,],tri)$rv)
    }
  }

  row.names(tri)<-c("vertex 1","vertex 2","vertex 3")  #vertex labeling

  list(rv=ind.set, #relative vertices
       tri=tri #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title Centers for non-degenerate asymptotic distribution of domination number of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#'
#' @description Returns the centers which yield nondegenerate asymptotic distribution for the domination number
#' of PE-PCD for uniform data in a triangle, \code{tri}\eqn{=T(v_1,v_2,v_3)}.
#'
#' PE proximity region is defined with
#' respect to the triangle \code{tri} with expansion parameter \code{r} in \eqn{(1,1.5]}.
#'
#' Vertex regions are defined with the centers that are output of this function.
#' Centers are stacked row-wise with row number is corresponding to the vertex
#' row number in \code{tri} (see the examples for an illustration). The center labels 1,2,3 correspond
#' to the vertices \eqn{M_1}, \eqn{M_2}, and \eqn{M_3} (which are the three centers for \code{r} in \eqn{(1,1.5)} which becomes center of mass
#' for \eqn{r=1.5}.).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]} for this function.
#'
#' @return The centers (stacked row-wise) which give nondegenerate asymptotic distribution
#' for the domination number of PE-PCD for uniform data in a triangle, \code{tri}.
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' r<-1.35
#'
#' Ms<-cent.nondegPE(Tr,r)
#' Ms
#'
#' Xlim<-range(Tr[,1])
#' Ylim<-range(Tr[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",main="Centers of nondegeneracy\n for the PE-PCD in a triangle",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Ms,pch=".",col=1)
#' polygon(Ms,lty=2)
#'
#' xc<-Tr[,1]+c(-.02,.02,.02)
#' yc<-Tr[,2]+c(.02,.02,.03)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' xc<-Ms[,1]+c(-.04,.04,.03)
#' yc<-Ms[,2]+c(.02,.02,.05)
#' txt.str<-c("M1","M2","M3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cent.nondegPE
cent.nondegPE <- function(tri,r)
{
  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!is.point(r,1) || r<=1 || r>1.5)
  {stop('r must be a scalar in (1,1.5]')}

  A<-tri[1,]; a1<-A[1]; a2<-A[2];
  B<-tri[2,]; b1<-B[1]; b2<-B[2];
  C<-tri[3,]; c1<-C[1]; c2<-C[2];

  M1x<--(a1*r-b1*r-c1*r-2*a1+b1+c1)/r;
  M1y<--(a2*r-b2*r-c2*r-2*a2+b2+c2)/r;
  M1<-c(M1x,M1y);

  M2x<-(a1*r-b1*r+c1*r-a1+2*b1-c1)/r;
  M2y<-(a2*r-b2*r+c2*r-a2+2*b2-c2)/r;
  M2<-c(M2x,M2y);

  M3x<-(a1*r+b1*r-c1*r-a1-b1+2*c1)/r;
  M3y<-(a2*r+b2*r-c2*r-a2-b2+2*c2)/r;
  M3<-c(M3x,M3y);

  rbind(M1,M2,M3)
} #end of the function
#'


#################################################################
#Auxiliary Function for Intervals
#################################################################

#' @title Parameterized center of an interval
#'
#' @description Returns the (parameterized) center, \eqn{M_c}, of the interval, \code{int}\eqn{=(a,b)}, parameterized by \eqn{c \in (0,1)}
#' so that \eqn{100c} \% of the length of interval is to the left of \eqn{M_c} and \eqn{100(1-c)} \% of the length of the interval
#' is to the right of \eqn{M_c}. That is, for the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param int A \code{vector} with two entries representing an interval.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return (parameterized) center inside \code{int}
#'
#' @seealso \code{\link{centersMc}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' c<-.4
#' a<-0; b<-10
#' int = c(a,b)
#' centMc(int,c)
#'
#' c<-.3
#' a<-2; b<-4; int<-c(a,b)
#' centMc(int,c)
#'
#' @export centMc
centMc <- function(int,c=.5)
{
  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  a<-int[1]; b<-int[2]
  Mc<-a+c*(b-a)
  Mc
} #end of the function
#'

#################################################################

#' @title Parameterized centers of intervals
#'
#' @description Returns the centers of the intervals based on 1D points in `Yp` parameterized by \eqn{c \in (0,1)}
#' so that \eqn{100c} \% of the length of interval is to the left of \eqn{M_c} and \eqn{100(1-c)} \% of the length of the interval
#' is to the right of \eqn{M_c}. That is, for an interval \eqn{(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}
#' \code{Yp} is a \code{vector} of 1D points, not necessarily sorted.
#'
#' See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param Yp A \code{vector} real numbers that constitute the end points of intervals.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the centers inside the intervals
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return (parameterized) centers of the intervals based on \code{Yp} points as a vector
#'
#' @seealso \code{\link{centMc}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-10
#' c<-.4  #try also c<-runif(1)
#' Yp<-runif(n)
#' centersMc(Yp,c)
#'
#' c<-.3  #try also c<-runif(1)
#' Yp<-runif(n,0,10)
#' centersMc(Yp,c)
#' }
#'
#' @export centersMc
centersMc <- function(Yp,c=.5)
{
  if (!is.point(Yp,length(Yp)))
  {stop('Yp must be a 1D vector of numerical entries')}

  if (!is.point(c,1) || c <= 0 || c >= 1 )
  {stop('c must be a scalar in (0,1)')}

  Yp<-sort(Yp)
  n<-length(Yp)
  Mcvec<-vector()
  for  (i in 1:(n-1))
  {int<-c(Yp[i],Yp[i+1])
  Mcvec<-c(Mcvec,centMc(int,c))
  }
  Mcvec
} #end of the function
#'

#' @title The index of the vertex region in a middle interval that contains a given point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' the interval \code{int}\eqn{=(a,b)=}(vertex 1,vertex 2) with (parameterized) center \eqn{M_c} associated with
#' the centrality parameter \eqn{c \in (0,1)}; vertices of interval are labeled as 1 and 2 according to their
#' order in the interval \code{int}. If the point, \code{p}, is not inside \code{int}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the interval \eqn{(a,b)} as \eqn{(a,M_c)} and \eqn{(M_c,b)} where \eqn{M_c=a+c(b-a)}.
#'
#' See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param p A 1D point. The vertex region \code{p} resides is to be found.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}
#' with the default \code{c=.5}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the vertex in the interval \code{int} whose region contains point, \code{p}.}
#' \item{int}{The vertices of the interval as a \code{vector} where position of the vertex corresponds to
#' the vertex index as \code{int=(rv=1,rv=2)}.}
#'
#' @seealso \code{\link{rv.end.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' a<-0; b<-10; int = c(a,b)
#'
#' Mc<-centMc(int,c)
#'
#' rv.mid.int(6,int,c)
#'
#' n<-20  #try also n<-40
#' xr<-range(a,b,Mc)
#' xf<-(int[2]-int[1])*.5
#' Xp<-runif(n,a,b)
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv,rv.mid.int(Xp[i],int,c)$rv)
#' Rv
#'
#' jit<-.1
#' yjit<-runif(n,-jit,jit)
#'
#' Xlim<-range(a,b,Xp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(Mc,0),main="vertex region indices for the points", xlab=" ", ylab=" ",
#'      xlim=Xlim+xd*c(-.05,.05),ylim=3*range(yjit),pch=".",cex=3)
#' abline(h=0)
#' points(Xp,yjit)
#' abline(v=c(a,b,Mc),lty=2,col=c(1,1,2))
#' text(Xp,yjit,labels=factor(Rv))
#' text(cbind(c(a,b,Mc),.02),c("rv=1","rv=2","Mc"))
#' }
#'
#' @export rv.mid.int
rv.mid.int <- function(p,int,c=.5)
{
  if (!is.point(p,1))
  {stop('p must be a scalar')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (p < y1 || p > y2)
  {rv<-NA
  } else
  {
    Mc<-y1+c*(y2-y1)
    rv <- ifelse(p<=Mc,1,2)
  }
  names(int)<-c("vertex 1","vertex 2")  #vertex labeling

  list(rv=rv, #relative vertex
       int=int #vertex labeling
  )
} #end of the function
#'

#' @title The index of the vertex region in an end-interval that contains a given point
#'
#' @description Returns the index of the vertex in the interval, \code{int}, whose end interval contains the 1D point \code{p},
#' that is, it finds the index of the vertex for the point, \code{p}, outside
#' the interval \code{int}\eqn{=(a,b)=}(vertex 1,vertex 2);
#' vertices of interval are labeled as 1 and 2 according to their order in the interval.
#'
#' If the point, \code{p}, is inside \code{int}, then the function yields \code{NA} as output.
#' The corresponding vertex region is an interval as \eqn{(-\infty,a)} or \eqn{(b,\infty)} for the interval \eqn{(a,b)}.
#' Then if \eqn{p<a}, then \code{rv=1} and if \eqn{p>b}, then \code{rv=2}. Unlike \code{\link{rv.mid.int}}, centrality parameter (i.e., center
#' of the interval is not relevant for \code{rv.end.int}.)
#'
#' See also (\insertCite{ceyhan:metrika-2012,ceyhan:revstat-2016;textual}{pcds}).
#'
#' @param p A 1D point whose end interval region is provided by the function.
#' @param int A \code{vector} of two real numbers representing an interval.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the end vertex whose region contains point, \code{p}.}
#' \item{int}{The vertices of the interval as a \code{vector} where position of the vertex corresponds to
#' the vertex index as \code{int=(rv=1,rv=2)}.}
#'
#' @seealso \code{\link{rv.mid.int}}
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
#' rv.end.int(-6,int)
#' rv.end.int(16,int)
#'
#' n<-10
#' xf<-(int[2]-int[1])*.5
#' XpL<-runif(n,a-xf,a)
#' XpR<-runif(n,b,b+xf)
#' Xp<-c(XpL,XpR)
#' rv.end.int(Xp[1],int)
#'
#' Rv<-vector()
#' for (i in 1:length(Xp))
#'   Rv<-c(Rv,rv.end.int(Xp[i],int)$rv)
#' Rv
#'
#' Xlim<-range(a,b,Xp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),xlab="",pch=".",xlim=Xlim+xd*c(-.05,.05))
#' abline(h=0)
#' abline(v=c(a,b),col=1,lty=2)
#' points(cbind(Xp,0))
#' text(cbind(Xp,0.1),labels=factor(Rv))
#' text(cbind(c(a,b),-0.1),c("rv=1","rv=2"))
#'
#' jit<-.1
#' yjit<-runif(length(Xp),-jit,jit)
#'
#' Xlim<-range(a,b,Xp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),main="vertex region indices for the points\n in the end intervals",
#'      xlab=" ", ylab=" ",pch=".",xlim=Xlim+xd*c(-.05,.05),ylim=3*range(yjit))
#' points(Xp, yjit,xlim=Xlim+xd*c(-.05,.05),pch=".",cex=3)
#' abline(h=0)
#' abline(v=c(a,b),lty=2)
#' text(Xp,yjit,labels=factor(Rv))
#' text(cbind(c(a,b),-.01),c("rv=1","rv=2"))
#' }
#'
#' @export rv.end.int
rv.end.int <- function(p,int)
{
  if (!is.point(p,1))
  {stop('p must be a scalar')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  if (p >= y1 & p <= y2)
  {stop('point must be outside the interval')}

  rv<-ifelse(p < y1,1,2)

  names(int)<-c("vertex 1","vertex 2")  #vertex labeling

  list(rv=rv, #relative vertex
       int=int #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The plot of the subintervals based on \code{Yp} points together with \code{Xp} points
#'
#' @description Plots the \code{Xp} points and the intervals based on \code{Yp} points
#' points.
#'
#' @param Xp A set of 1D points whose scatter-plot is provided.
#' @param Yp A set of 1D points which constitute the end points of the intervals which
#' partition the real line.
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the intervals based on \code{Yp} points and also scatter plot of \code{Xp} points
#'
#' @seealso \code{\link{plotPEregs1D}} and \code{\link{plotDeltri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' a<-0; b<-10;
#'
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-4;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-runif(nx,a,b)
#' Yp<-runif(ny,a,b)
#'
#' plotIntervals(Xp,Yp,xlab="",ylab="")
#' }
#'
#' @export plotIntervals
plotIntervals <- function(Xp,Yp,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL, ...)
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.point(Xp,length(Xp)) || !is.point(Yp,length(Yp)) )
  {stop('Xp and Yp must be 1D vectors of numerical entries')}

  nx<-length(Xp); ny<-length(Yp)

  if (ny<1 || nx<1)
  {stop('Both Xp and Yp points must be nonempty to construct intervals')}

  Ys<-sort(Yp)
  LE<-Ys[1:(ny-1)]
  RE<-Ys[2:ny]

  if (is.null(xlim))
  {xl<-range(Xp,Yp)
  xr<-xl[2]-xl[1]
  xlim<-xl+.05*xr*c(-1,1)
  } else
  {
    xr<-xlim[2]-xlim[1]
    xlim<-xlim+.05*xr*c(-1,1)
  }

  if (is.null(ylim))
  {yl<-xr*.05*sin(30*pi/180)
  ylim<-yl*c(-1,1)
  }

  if (is.null(main))
  { main=paste(" Points in ", xname," and Intervals based on Points in ",yname,sep="")}

  plot(cbind(Xp, 0),main=main, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3, ...)
  points(cbind(Yp,0), col=2, ...)
  abline(h=0,lty=1)
  for (i in 1:ny)
  {
    plotrix::draw.arc(LE[i]+xr*.1, 0,xr*.1, deg1=150,deg2 = 210, col = "blue")
    plotrix::draw.arc(RE[i]-xr*.1, 0, xr*.1, deg1=-30,deg2 = 30, col = "blue")
  }
} #end of the function
#'


#################################################################
#Auxiliary Function for Tetrahedrons
#################################################################

#' @title Check whether a point is inside a tetrahedron
#'
#' @description Checks if the point \code{p} lies in the tetrahedron, \code{th}, using the barycentric coordinates, generally denoted as
#' \eqn{(\alpha,\beta,\gamma)}. If all (normalized or non-normalized) barycentric coordinates are positive then the
#' point \code{p} is inside the tetrahedron, if all are nonnegative with one or more are zero, then \code{p} falls on the
#' boundary. If some of the barycentric coordinates are negative, then \code{p} falls outside the tetrahedron.
#'
#' \code{boundary} is a logical argument (default=\code{FALSE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if the point, \code{p}, lies in the closure of the tetrahedron (i.e., interior and boundary
#' combined) else it checks if \code{p} lies in the interior of the tetrahedron.
#'
#' @param p A 3D point to be checked whether it is inside the tetrahedron or not.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param boundary A logical parameter (default=\code{TRUE}) to include boundary or not, so if it is \code{TRUE},
#' the function checks if the point, \code{p}, lies in the closure of the tetrahedron (i.e., interior and boundary
#' combined); else, it checks if \code{p} lies in the interior of the tetrahedron.
#'
#' @return A \code{list} with two elements
#' \item{in.tetra}{A logical output, if the point, \code{p}, is inside the tetrahedron, \code{th}, it is \code{TRUE},
#' else it is \code{FALSE}.}
#' \item{barycentric}{The barycentric coordinates of the point \code{p} with respect to the tetrahedron, \code{th}.}
#'
#' @seealso \code{\link{in.triangle}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3); P<-c(.1,.1,.1)
#' tetra<-rbind(A,B,C,D)
#'
#' in.tetrahedron(P,tetra,boundary=FALSE)
#'
#' in.tetrahedron(C,tetra)
#' in.tetrahedron(C,tetra,boundary=FALSE)
#'
#' n1<-5; n2<-5; n<-n1+n2
#' Xp<-rbind(cbind(runif(n1),runif(n1,0,sqrt(3)/2),runif(n1,0,sqrt(6)/3)),
#'           runif.tetra(n2,tetra)$g)
#'
#' in.tetra<-vector()
#' for (i in 1:n)
#' {in.tetra<-c(in.tetra,in.tetrahedron(Xp[i,],tetra,boundary=TRUE)$in.tetra) }
#'
#' in.tetra
#' dat.tet<-Xp[in.tetra,]
#' if (is.vector(dat.tet)) {dat.tet<-matrix(dat.tet,nrow=1)}
#'
#' Xlim<-range(tetra[,1],Xp[,1])
#' Ylim<-range(tetra[,2],Xp[,2])
#' Zlim<-range(tetra[,3],Xp[,3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi=40,theta=40, bty = "g", pch = 20, cex = 1,
#' ticktype="detailed",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.05,.05))
#' #add the vertices of the tetrahedron
#' plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
#' plot3D::points3D(dat.tet[,1],dat.tet[,2],dat.tet[,3],pch=4, add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#'
#' in.tetrahedron(P,tetra) #this works fine
#' }
#'
#' @export in.tetrahedron
in.tetrahedron <- function(p,th,boundary=TRUE)
{
  if (!is.point(p,3))
  {stop('p must be a numeric 3D point')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (isTRUE(all.equal(D0, 0)) )
  {stop('the tetrahedron is degenerate')}

  if (p[1]=="NA" || p[2]=="NA"|| p[3]=="NA")
  {ind.tetra<-FALSE; b1<-b2<-b3<-b4<-NA;
  } else
  {
    V1<-th[1,];V2<-th[2,];V3<-th[3,];V4<-th[4,];

    D1<-det(matrix(cbind(rbind(p,V2,V3,V4),vec1),ncol=4))
    D2<-det(matrix(cbind(rbind(V1,p,V3,V4),vec1),ncol=4))
    D3<-det(matrix(cbind(rbind(V1,V2,p,V4),vec1),ncol=4))
    D4<-det(matrix(cbind(rbind(V1,V2,V3,p),vec1),ncol=4))
    b1<-D1/D0; b2<-D2/D0; b3<-D3/D0; b4<-D4/D0; # barycentric coordinates

    b.c = c(b1,b2,b3,b4)
    ind0 = abs(b.c) <= .Machine$double.eps #indices of the zeros
    b.c[ind0] = 0

    #if all bi are greater than 0, then the point p lies within the tetrahedron
    if (boundary==TRUE)
    {ifelse(all(b.c>=0), ind.tetra<-TRUE, ind.tetra<-FALSE)
    }  else
    {ifelse(all(b.c>0), ind.tetra<-TRUE, ind.tetra<-FALSE)}
  }

  list(
    in.tetra=ind.tetra,
    barycentric=b.c
  )
}
#'

#################################################################


#' @title Circumcenter of a general tetrahedron
#'
#' @description Returns the circumcenter a given tetrahedron \code{th} with vertices stacked row-wise.
#'
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#'
#' @return circumcenter of the tetrahedron \code{th}
#'
#' @seealso \code{\link{circ.cent.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' A<-c(0,0,0)+runif(3,-.2,.2);
#' B<-c(1,0,0)+runif(3,-.2,.2);
#' C<-c(1/2,sqrt(3)/2,0)+runif(3,-.2,.2);
#' D<-c(1/2,sqrt(3)/6,sqrt(6)/3)+runif(3,-.2,.2);
#' tetra<-rbind(A,B,C,D)
#'
#' CC<-circ.cent.tetra(tetra)
#' CC
#'
#' Xlim<-range(tetra[,1],CC[1])
#' Ylim<-range(tetra[,2],CC[2])
#' Zlim<-range(tetra[,3],CC[3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(tetra[,1],tetra[,2],tetra[,3], phi =0,theta=40, bty = "g",
#' main="Illustration of the Circumcenter\n in a Tetrahedron",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
#'           pch = 20, cex = 1, ticktype = "detailed")
#' #add the vertices of the tetrahedron
#' plot3D::points3D(CC[1],CC[2],CC[3], add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#'
#' D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
#' L<-rbind(D1,D2,D3,D4,D5,D6); R<-matrix(rep(CC,6),byrow = TRUE,ncol=3)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
#'
#' plot3D::text3D(CC[1],CC[2],CC[3], labels="CC", add=TRUE)
#' }
#'
#' @export circ.cent.tetra
circ.cent.tetra <- function(th)
{
  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('the argument must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  A<-th[1,]; B<-th[2,]; C<-th[3,]; D<-th[4,];
  a1<-A[1]; a2<-A[2]; a3<-A[3];
  b1<-B[1]; b2<-B[2]; b3<-B[3];
  c1<-C[1]; c2<-C[2]; c3<-C[3];

  Dx<-det(matrix(cbind(c(sum(A^2),sum(B^2),sum(C^2),sum(D^2)),th[,2],th[,3],vec1),ncol=4))
  Dy<--det(matrix(cbind(c(sum(A^2),sum(B^2),sum(C^2),sum(D^2)),th[,1],th[,3],vec1),ncol=4))
  Dz<-det(matrix(cbind(c(sum(A^2),sum(B^2),sum(C^2),sum(D^2)),th[,1],th[,2],vec1),ncol=4))

  cc1<-Dx/(2*D0); cc2<-Dy/(2*D0); cc3<-Dz/(2*D0);
  c(cc1,cc2,cc3)
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CM}-vertex region in a tetrahedron that contains a point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' a tetrahedron \eqn{th=T(A,B,C,D)} and vertex regions are based on the center of mass \eqn{CM=(A+B+C+D)/4} of \code{th}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the tetrahedron \code{th} are labeled as
#' \eqn{1=A}, \eqn{2=B}, \eqn{3=C}, and \eqn{4=C} also according to the row number the vertex is recorded in \code{th}.
#'
#' If the point, \code{p}, is not inside \code{th}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the simplex with the vertex, \eqn{CM}, and
#' midpoints of the edges adjacent to the vertex.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p A 3D point for which \eqn{CM}-vertex region it resides in is to be determined in the
#' tetrahedron \code{th}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the \eqn{CM}-vertex region that contains point, \code{p} in the tetrahedron \code{th}}
#' \item{th}{The vertices of the tetrahedron, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rv.tetraCC}} and \code{\link{rv.triCM}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' tetra<-rbind(A,B,C,D)
#'
#' n<-20  #try also n<-40
#'
#' Xp<-runif.std.tetra(n)$g
#'
#' rv.tetraCM(Xp[1,],tetra)
#'
#' Rv<-vector()
#' for (i in 1:n)
#'   Rv<-c(Rv, rv.tetraCM(Xp[i,],tetra)$rv )
#' Rv
#'
#' Xlim<-range(tetra[,1],Xp[,1])
#' Ylim<-range(tetra[,2],Xp[,2])
#' Zlim<-range(tetra[,3],Xp[,3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' CM<-apply(tetra,2,mean)
#'
#' plot3D::scatter3D(tetra[,1],tetra[,2],tetra[,3], phi =0,theta=40, bty = "g",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
#'           pch = 20, cex = 1, ticktype = "detailed")
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#' #add the data points
#' plot3D::points3D(Xp[,1],Xp[,2],Xp[,3],pch=".",cex=3, add=TRUE)
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#' plot3D::text3D(CM[1],CM[2],CM[3], labels=c("CM"), add=TRUE)
#'
#' D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
#' L<-rbind(D1,D2,D3,D4,D5,D6); R<-matrix(rep(CM,6),ncol=3,byrow=TRUE)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
#'
#' F1<-int.line.plane(A,CM,B,C,D)
#' L<-matrix(rep(F1,4),ncol=3,byrow=TRUE); R<-rbind(D4,D5,D6,CM)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=2, add=TRUE,lty=2)
#'
#' F2<-int.line.plane(B,CM,A,C,D)
#' L<-matrix(rep(F2,4),ncol=3,byrow=TRUE); R<-rbind(D2,D3,D6,CM)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=3, add=TRUE,lty=2)
#'
#' F3<-int.line.plane(C,CM,A,B,D)
#' L<-matrix(rep(F3,4),ncol=3,byrow=TRUE); R<-rbind(D3,D5,D6,CM)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=4, add=TRUE,lty=2)
#'
#' F4<-int.line.plane(D,CM,A,B,C)
#' L<-matrix(rep(F4,4),ncol=3,byrow=TRUE); R<-rbind(D1,D2,D4,CM)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=5, add=TRUE,lty=2)
#'
#' plot3D::text3D(Xp[,1],Xp[,2],Xp[,3], labels=factor(Rv), add=TRUE)
#' }
#'
#' @export rv.tetraCM
rv.tetraCM <- function(p,th)
{
  if (!is.point(p,3))
  {stop('p must be a numeric 3D point')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (!in.tetrahedron(p,th,boundary=TRUE)$in.tetra)
  {rv<-NA
  } else
  {
    CM<-apply(th,2,mean)
    i<-1
    while (i<=4)
    { th0<-th[-i,] #remove ith vertex
    y1<-th[i,]; y2<-th0[1,]; y3<-th0[2,]; y4<-th0[3,];
    M1<-(y1+y2)/2; M2<-(y1+y3)/2; M3<-(y1+y4)/2;
    th1<-rbind(y1,M1,y3,y4); th2<-rbind(y1,M2,y2,y4); th3<-rbind(y1,M3,y2,y3);
    if (in.tetrahedron(p,th1,boundary=TRUE)$i &&
        in.tetrahedron(p,th2,boundary=TRUE)$i &&
        in.tetrahedron(p,th3,boundary=TRUE)$i)
    {rv<-i; i<-5
    } else
    {i<-i+1}
    }
  }
  row.names(th)<-c("vertex 1","vertex 2","vertex 3","vertex 4")  #vertex labeling

  list(rv=rv, #relative vertex
       tetra=th #vertex labeling
  )
} #end of the function
#'

#################################################################

#' @title The index of the \eqn{CC}-vertex region in a tetrahedron that contains a point
#'
#' @description Returns the index of the vertex whose region contains point \code{p} in
#' a tetrahedron \eqn{th=T(A,B,C,D)} and vertex regions are based on the circumcenter \eqn{CC} of \code{th}.
#' (see the plots in the example for illustrations).
#'
#' The vertices of the tetrahedron \code{th} are labeled as
#' \eqn{1=A}, \eqn{2=B}, \eqn{3=C}, and \eqn{4=C} also according to the row number the vertex is recorded in \code{th}.
#'
#' If the point, \code{p}, is not inside \code{th}, then the function yields \code{NA} as output.
#' The corresponding vertex region is the polygon whose interior points are closest to that vertex.
#' If \code{th} is regular tetrahedron, then \eqn{CC} and \eqn{CM} (center of mass) coincide.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p A 3D point for which \eqn{CC}-vertex region it resides in is to be determined in the
#' tetrahedron \code{th}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#'
#' @return A \code{list} with two elements
#' \item{rv}{Index of the \eqn{CC}-vertex region that contains point, \code{p} in the tetrahedron \code{th}}
#' \item{tri}{The vertices of the tetrahedron, where row number corresponds to the vertex index in \code{rv}.}
#'
#' @seealso \code{\link{rv.tetraCM}} and \code{\link{rv.triCC}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' A<-c(0,0,0)+runif(3,-.2,.2);
#' B<-c(1,0,0)+runif(3,-.2,.2);
#' C<-c(1/2,sqrt(3)/2,0)+runif(3,-.2,.2);
#' D<-c(1/2,sqrt(3)/6,sqrt(6)/3)+runif(3,-.2,.2);
#' tetra<-rbind(A,B,C,D)
#'
#' n<-20  #try also n<-40
#'
#' Xp<-runif.tetra(n,tetra)$g
#'
#' rv.tetraCC(Xp[1,],tetra)
#'
#' Rv<-vector()
#' for (i in 1:n)
#'  Rv<-c(Rv,rv.tetraCC(Xp[i,],tetra)$rv)
#' Rv
#'
#' CC<-circ.cent.tetra(tetra)
#' CC
#'
#' Xlim<-range(tetra[,1],Xp[,1],CC[1])
#' Ylim<-range(tetra[,2],Xp[,2],CC[2])
#' Zlim<-range(tetra[,3],Xp[,3],CC[3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(tetra[,1],tetra[,2],tetra[,3], phi =0,theta=40, bty = "g",
#' main="Scatterplot of data points \n and CC-vertex regions",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
#'           pch = 20, cex = 1, ticktype = "detailed")
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#' #add the data points
#' plot3D::points3D(Xp[,1],Xp[,2],Xp[,3],pch=".",cex=3, add=TRUE)
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#' plot3D::text3D(CC[1],CC[2],CC[3], labels=c("CC"), add=TRUE)
#'
#' D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
#' L<-rbind(D1,D2,D3,D4,D5,D6); R<-matrix(rep(CC,6),ncol=3,byrow=TRUE)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
#'
#' F1<-int.line.plane(A,CC,B,C,D)
#' L<-matrix(rep(F1,4),ncol=3,byrow=TRUE); R<-rbind(D4,D5,D6,CC)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=2, add=TRUE,lty=2)
#'
#' F2<-int.line.plane(B,CC,A,C,D)
#' L<-matrix(rep(F2,4),ncol=3,byrow=TRUE); R<-rbind(D2,D3,D6,CC)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=3, add=TRUE,lty=2)
#'
#' F3<-int.line.plane(C,CC,A,B,D)
#' L<-matrix(rep(F3,4),ncol=3,byrow=TRUE); R<-rbind(D3,D5,D6,CC)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=4, add=TRUE,lty=2)
#'
#' F4<-int.line.plane(D,CC,A,B,C)
#' L<-matrix(rep(F4,4),ncol=3,byrow=TRUE); R<-rbind(D1,D2,D4,CC)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3],col=5, add=TRUE,lty=2)
#'
#' plot3D::text3D(Xp[,1],Xp[,2],Xp[,3], labels=factor(Rv), add=TRUE)
#' }
#'
#' @export rv.tetraCC
rv.tetraCC <- function(p,th)
{
  if (!is.point(p,3))
  {stop('p must be a numeric 3D point')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (!in.tetrahedron(p,th,boundary=TRUE)$in.tetra)
  {rv<-NA
  } else
  {
    y1<-th[1,]; y2<-th[2,]; y3<-th[3,]; y4<-th[4,];
    mdt<-max(Dist(y1,y2),Dist(y1,y3),Dist(y1,y4),Dist(y2,y3),Dist(y2,y4),Dist(y3,y4))
    for (i in 1:4)
    {
      d1<-Dist(p,th[i,]);
      if (d1<mdt)
      {mdt<-d1; rv<-i }
    }
  }
  row.names(th)<-c("vertex 1","vertex 2","vertex 3","vertex 4")  #vertex labeling

  list(rv=rv, #relative vertex
       tetra=th #vertex labeling
  )
} #end of the function
#'
