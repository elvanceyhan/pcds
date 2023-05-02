#ArcSliceFunctions.R;
#Functions for NAS in R^2

#################################################################

# funsTbMid2CC
#'
#' @title Two functions \code{lD1CCinTb} and \code{lD2CCinTb} which are of class \code{"TriLines"} ---
#' The lines joining the midpoints of edges to the circumcenter (\eqn{CC}) in the standard basic triangle.
#'
#' @description Returns the \code{equation, slope, intercept}, and \eqn{y}-coordinates of the lines
#' joining \eqn{D_1} and \eqn{CC} and also \eqn{D_2} and \eqn{CC}, in the standard basic triangle \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}
#' and \eqn{D_1=(B+C)/2} and \eqn{D_2=(A+C)/2} are the midpoints of edges \eqn{BC} and \eqn{AC}.
#'
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#' \eqn{x}-coordinates are provided in \code{vector} \code{x}.
#'
#' @param x A single scalar or a \code{vector} of scalars.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Longer description of the line.}
#' \item{txt2}{Shorter description of the line (to be inserted over the line in the plot).}
#' \item{mtitle}{The \code{"main"} title for the plot of the line.}
#' \item{cent}{The center chosen inside the standard equilateral triangle.}
#' \item{cent.name}{The name of the center inside the standard basic triangle.
#' It is \code{"CC"} for these two functions.}
#' \item{tri}{The triangle (it is the standard basic triangle for this function).}
#' \item{x}{The input vector, can be a scalar or a \code{vector} of scalars,
#' which constitute the \eqn{x}-coordinates of the point(s) of interest on the line.}
#' \item{y}{The output vector, will be a scalar if \code{x} is a scalar or a \code{vector} of scalars if \code{x} is a \code{vector} of scalar,
#' constitutes the \eqn{y}-coordinates of the point(s) of interest on the line.}
#' \item{slope}{Slope of the line.}
#' \item{intercept}{Intercept of the line.}
#' \item{equation}{Equation of the line.}
#'
#' @name funsTbMid2CC
NULL
#'
#' @seealso \code{\link{lA_CM.Te}}, \code{\link{lB_CM.Te}}, \code{\link{lA_M.Te}}, \code{\link{lB_M.Te}},
#' and \code{\link{lC_M.Te}}
#'
#' @rdname funsTbMid2CC
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for lD1CCinTb
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);  #the vertices of the standard basic triangle Tb
#'
#' Tb<-rbind(A,B,C)
#'
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.1)  #try also by=.01
#'
#' lnD1CC<-lD1CCinTb(x,c1,c2)
#' lnD1CC
#' summary(lnD1CC)
#' plot(lnD1CC)
#'
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
#' Ds<-rbind(D1,D2,D3)
#'
#' x1<-seq(0,1,by=.1)  #try also by=.01
#' y1<-lD1CCinTb(x1,c1,c2)$y
#'
#' Xlim<-range(Tb[,1],x1)
#' Ylim<-range(Tb[,2],y1)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tb,CC,D1,D2,D3)
#' xc<-txt[,1]+c(-.03,.04,.03,.02,.09,-.08,0)
#' yc<-txt[,2]+c(.02,.02,.04,.08,.03,.03,-.05)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' lines(x1,y1,type="l",lty=2)
#' text(.8,.5,"lD1CCinTb")
#'
#' c1<-.4; c2<-.6;
#' x1<-seq(0,1,by=.1)  #try also by=.01
#' lD1CCinTb(x1,c1,c2)
#' }
#'
#' @export
lD1CCinTb <- function(x,c1,c2)
{
  dname <-deparse(substitute(x))

  if (!is.vector(x))
  {stop('x must be a vector')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  A<-c(0,0); B<-c(1,0); C<-c(c1,c2);  #the vertices of the standard basic triangle
  Tr<-rbind(A,B,C)
  M<-circ.cent.tri(Tr)

  sl<-(1-c1)/c2
  int<--1/2*(-c2^2-c1^2-2*0+2*c1*0+1)/c2
  ln<-int+sl*x

  names(sl)<-"slope"
  names(int)<-"intercept"
  cname <-"CC"

  txt1<-"Line joining midpoint of the edge BC, D1, to CC (circumcenter) in the standard basic triangle Tb=T(A,B,C)=T((0,0),(1,0),(c1,c2))"
  txt2<-paste("lD1CCinTb(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2,
    cent=M, cent.name=cname,
    tri=Tr,
    x=x,
    y=ln,
    slope=sl,
    intercept=int,
    equation=paste("y=",sl,"x+",int,sep="")
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'
#' @rdname funsTbMid2CC
#'
#' @examples
#' \dontrun{
#' #Examples for lD2CCinTb
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);  #the vertices of the standard basic triangle Tb
#'
#' Tb<-rbind(A,B,C)
#'
#' xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,by=.1)  #try also by=.01
#'
#' lnD2CC<-lD2CCinTb(x,c1,c2)
#' lnD2CC
#' summary(lnD2CC)
#' plot(lnD2CC)
#'
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
#' Ds<-rbind(D1,D2,D3)
#'
#' x2<-seq(0,1,by=.1)  #try also by=.01
#' y2<-lD2CCinTb(x2,c1,c2)$y
#'
#' Xlim<-range(Tb[,1],x1)
#' Ylim<-range(Tb[,2],y2)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' txt<-rbind(Tb,CC,D1,D2,D3)
#' xc<-txt[,1]+c(-.03,.04,.03,.02,.09,-.08,0)
#' yc<-txt[,2]+c(.02,.02,.04,.08,.03,.03,-.05)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' lines(x2,y2,type="l",lty=2)
#' text(0,.5,"lD2CCinTb")
#' }
#'
#' @export
lD2CCinTb <- function(x,c1,c2)
{
  dname <-deparse(substitute(x))

  if (!is.vector(x))
  {stop('x must be a vector')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  A<-c(0,0); B<-c(1,0); C<-c(c1,c2);  #the vertices of the standard basic triangle
  Tr<-rbind(A,B,C)
  M<-circ.cent.tri(Tr)

  sl<--c1/c2
  int<--1/2*(-c2^2-c1^2+2*c1*0)/c2
  ln<-int+sl*x

  names(sl)<-"slope"
  names(int)<-"intercept"
  cname <-"CC"

  txt1<-"Line joining midpoint of the edge AC, D2, to CC (circumcenter) in the standard basic triangle Tb=T(A,B,C)=T((0,0),(1,0),(c1,c2))"
  txt2<-paste("lD2CCinTb(",dname,")",sep="")

  res<-list(
    txt1=txt1, txt2=txt2,
    cent=M, cent.name=cname,
    tri=Tr,
    x=x,
    y=ln,
    slope=sl,
    intercept=int,
    equation=paste("y=",sl,"x+",int,sep="")
  )

  class(res)<-"TriLines"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The angle between two line segments
#'
#' @description Returns the angle in radians or degrees between two vectors or line segments at the point of
#' intersection. The line segments are \eqn{[ba]} and \eqn{[bc]} when the arguments of the function are given as \code{a,b,c}.
#' \code{radian} is a logical argument (default=\code{TRUE}) which yields the angle in radians if \code{TRUE}, and in degrees if \code{FALSE}.
#' The smaller of the angle between the line segments is provided by the function.
#'
#' @param a,b,c Three 2D points which represent the intersecting line segments \eqn{[ba]} and \eqn{[bc]}.
#' The smaller angle between these line segments is to be computed.
#' @param radian A logical argument (default=\code{TRUE}). If \code{TRUE}, the (smaller) angle between the line segments
#' \eqn{[ba]} and \eqn{[bc]} is provided in radians, else it is provided in degrees.
#'
#' @return angle in radians or degrees between the line segments \eqn{[ba]} and \eqn{[bc]}
#'
#' @seealso \code{\link{angle.str2end}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(.3,.2); B<-c(.6,.3); C<-c(1,1)
#' pts<-rbind(A,B,C)
#'
#' angle3pnts(A,B,C)
#'
#' angle3pnts(A,B,A)
#' round(angle3pnts(A,B,A),7)
#'
#' angle3pnts(A,B,C,radian=FALSE)
#'
#' #plot of the line segments
#' Xlim<-range(pts[,1])
#' Ylim<-range(pts[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' ang1<-angle3pnts(A,B,C,radian=FALSE)
#' ang2<-angle3pnts(B+c(1,0),B,C,radian=FALSE)
#'
#' sa<-angle.str2end(A,B,C,radian=FALSE)$s #small arc angles
#' ang1<-sa[1]
#' ang2<-sa[2]
#'
#' plot(pts,asp=1,pch=1,xlab="x",ylab="y",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' L<-rbind(B,B); R<-rbind(A,C)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' plotrix::draw.arc(B[1],B[2],radius=xd*.1,deg1=ang1,deg2=ang2)
#' txt<-rbind(A,B,C)
#' text(txt+cbind(rep(xd*.05,nrow(txt)),rep(-xd*.02,nrow(txt))),c("A","B","C"))
#'
#' text(rbind(B)+.15*xd*c(cos(pi*(ang2+ang1)/360),sin(pi*(ang2+ang1)/360)),
#' paste(round(abs(ang1-ang2),2)," degrees"))
#' }
#'
#' @export angle3pnts
angle3pnts <- function(a,b,c,radian=TRUE)
{
  if (!is.point(a) || !is.point(b) || !is.point(c) )
  {stop('all entries must be all numeric 2D points')}

  u<-a-b; v<-c-b
  rhs<-sum(u*v)/(sqrt(sum(u*u))*sqrt(sum(v*v)))

  if  (radian==T)
  {angle<-acos(rhs)  #in radians
  } else
  {
    angle<-acos(rhs)*180/pi # in degrees
  }
  angle
} #end of the function
#'

#################################################################

#' @title The angles to draw arcs between two line segments
#'
#' @description
#' Gives the two pairs of angles in radians or degrees to draw arcs between two vectors or line segments
#' for the \code{\link[plotrix]{draw.arc}} function in the \code{plotrix} package.
#' The angles are provided with respect to the \eqn{x}-axis in the coordinate system.
#' The line segments are \eqn{[ba]} and \eqn{[bc]} when the argument is given as \code{a,b,c} in the function.
#'
#' \code{radian} is a logical argument (default=\code{TRUE}) which yields the angle in radians if \code{TRUE}, and in degrees if \code{FALSE}.
#' The first pair of angles is for drawing arcs in the smaller angle between \eqn{[ba]} and \eqn{[bc]}
#' and the second pair of angles is for drawing arcs in the counter-clockwise order from \eqn{[ba]} to \eqn{[bc]}.
#'
#' @param a,b,c Three 2D points which represent the intersecting line segments \eqn{[ba]} and \eqn{[bc]}.
#' @param radian A logical argument (default=\code{TRUE}).
#' If \code{TRUE}, the smaller angle or counter-clockwise angle
#' between the line segments \eqn{[ba]} and \eqn{[bc]} is provided in radians, else it is provided in degrees.
#'
#' @return A \code{list} with two elements
#' \item{small.arc.angles}{Angles of \eqn{[ba]} and \eqn{[bc]} with the \eqn{x}-axis so that difference between them
#' is the smaller angle between \eqn{[ba]} and \eqn{[bc]} }
#' \item{ccw.arc.angles}{Angles of \eqn{[ba]} and \eqn{[bc]} with the \eqn{x}-axis so that difference between them
#' is the counter-clockwise angle between \eqn{[ba]} and \eqn{[bc]}}
#'
#' @seealso \code{\link{angle3pnts}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(.3,.2); B<-c(.6,.3); C<-c(1,1)
#'
#' pts<-rbind(A,B,C)
#'
#' Xp<-c(B[1]+max(abs(C[1]-B[1]),abs(A[1]-B[1])),0)
#'
#' angle.str2end(A,B,C)
#' angle.str2end(A,B,A)
#'
#' angle.str2end(A,B,C,radian=FALSE)
#'
#' #plot of the line segments
#' ang.rad<-angle.str2end(A,B,C,radian=TRUE); ang.rad
#' ang.deg<-angle.str2end(A,B,C,radian=FALSE); ang.deg
#' ang.deg1<-ang.deg$s; ang.deg1
#' ang.deg2<-ang.deg$c; ang.deg2
#'
#' rad<-min(Dist(A,B),Dist(B,C))
#'
#' Xlim<-range(pts[,1],Xp[1],B+Xp,B[1]+c(+rad,-rad))
#' Ylim<-range(pts[,2],B[2]+c(+rad,-rad))
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' #plot for the smaller arc
#' plot(pts,pch=1,asp=1,xlab="x",ylab="y",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' L<-rbind(B,B,B); R<-rbind(A,C,B+Xp)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' plotrix::draw.arc(B[1],B[2],radius=.3*rad,angle1=ang.rad$s[1],angle2=ang.rad$s[2])
#' plotrix::draw.arc(B[1],B[2],radius=.6*rad,angle1=0, angle2=ang.rad$s[1],lty=2,col=2)
#' plotrix::draw.arc(B[1],B[2],radius=.9*rad,angle1=0,angle2=ang.rad$s[2],col=3)
#' txt<-rbind(A,B,C)
#' text(txt+cbind(rep(xd*.02,nrow(txt)),rep(-xd*.02,nrow(txt))),c("A","B","C"))
#'
#' text(rbind(B)+.5*rad*c(cos(mean(ang.rad$s)),sin(mean(ang.rad$s))),
#'      paste(abs(round(ang.deg1[2]-ang.deg1[1],2))," degrees",sep=""))
#' text(rbind(B)+.6*rad*c(cos(ang.rad$s[1]/2),sin(ang.rad$s[1]/2)),paste(round(ang.deg1[1],2)),col=2)
#' text(rbind(B)+.9*rad*c(cos(ang.rad$s[2]/2),sin(ang.rad$s[2]/2)),paste(round(ang.deg1[2],2)),col=3)
#'
#' #plot for the counter-clockwise arc
#' plot(pts,pch=1,asp=1,xlab="x",ylab="y",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' L<-rbind(B,B,B); R<-rbind(A,C,B+Xp)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' plotrix::draw.arc(B[1],B[2],radius=.3*rad,angle1=ang.rad$c[1],angle2=ang.rad$c[2])
#' plotrix::draw.arc(B[1],B[2],radius=.6*rad,angle1=0, angle2=ang.rad$s[1],lty=2,col=2)
#' plotrix::draw.arc(B[1],B[2],radius=.9*rad,angle1=0,angle2=ang.rad$s[2],col=3)
#' txt<-pts
#' text(txt+cbind(rep(xd*.02,nrow(txt)),rep(-xd*.02,nrow(txt))),c("A","B","C"))
#'
#' text(rbind(B)+.5*rad*c(cos(mean(ang.rad$c)),sin(mean(ang.rad$c))),
#'      paste(abs(round(ang.deg2[2]-ang.deg2[1],2))," degrees",sep=""))
#' text(rbind(B)+.6*rad*c(cos(ang.rad$s[1]/2),sin(ang.rad$s[1]/2)),paste(round(ang.deg1[1],2)),col=2)
#' text(rbind(B)+.9*rad*c(cos(ang.rad$s[2]/2),sin(ang.rad$s[2]/2)),paste(round(ang.deg1[2],2)),col=3)
#' }
#'
#' @export angle.str2end
angle.str2end <- function(a,b,c,radian=TRUE)
{
  if (!is.point(a) || !is.point(b) || !is.point(c) )
  {stop('all entries must be all numeric 2D points')}

  if (isTRUE(all.equal(a,b))==T || isTRUE(all.equal(c,b))==T)
  {stop('The angle is undefined')}

  u<-a-b; v<-c-b
  a1<-atan2(u[2],u[1])
  a2<-atan2(v[2],v[1])

  if (sign(a1)!=sign(a2))
  {
    if (a2<0) {A1<-a1; A2<-2*pi+a2} else {A1<-a1; A2<-a2}
    if (abs(a1)+abs(a2)<pi)
    { st.ang<-min(a1,a2);
    end.ang<-max(a1,a2)
    } else {
      B1<-max(a1,a2)
      B2<-2*pi+min(a1,a2)
      st.ang<-min(B1,B2);
      end.ang<-max(B1,B2)
    }
  } else
  {
    if (a1>a2) {A1<-a1; A2<-a2+2*pi} else {A1<-a1; A2<-a2}
    st.ang<-min(a1,a2);
    end.ang<-max(a1,a2)
  }

  if  (radian==T)
  {angles<-c(st.ang,end.ang)  #in radians
  } else
  {
    angles<-c(st.ang,end.ang)*180/pi # in degrees
    A1<-A1*180/pi; A2<-A2*180/pi;
  }
  list(
    small.arc.angles=angles, #angles are given so that arc between [ba] and [bc] is plotted in the smaller angle
    ccw.arc.angles=c(A1,A2)  #angles are given so that arc between [ba] and [bc] is plotted in counter-clockwise order,
  )
} #end of the function
#'

#################################################################

#' @title The points of intersection of a line and a circle
#'
#' @description Returns the intersection point(s) of a line and a circle. The line is determined by the two points
#' \code{p1} and \code{p2} and the circle is centered at point \code{cent} and has radius \code{rad}.
#' If the circle does not intersect the line, the function yields \code{NULL};
#' if the circle intersects at only one point, it yields only that point;
#' otherwise it yields both intersection points as output.
#' When there are two intersection points, they are listed in the order of the \eqn{x}-coordinates of \code{p1} and \code{p2};
#' and if the \eqn{x}-coordinates of \code{p1} and \code{p2} are equal, intersection points are listed in the order of
#' \eqn{y}-coordinates of \code{p1} and \code{p2}.
#'
#' @param p1,p2 2D points that determine the straight line (i.e., through which the straight line passes).
#' @param cent A 2D point representing the center of the circle.
#' @param rad A positive real number representing the radius of the circle.
#'
#' @return point(s) of intersection between the circle and the line (if they do not intersect, the function
#' yields \code{NULL} as the output)
#'
#' @seealso \code{\link{int.2lines}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' P1<-c(.3,.2)*100
#' P2<-c(.6,.3)*100
#' cent<-c(1.1,1.1)*100
#' rad<-2*100
#'
#' int.circ.line(P1,P2,cent,rad)
#' int.circ.line(P2,P1,cent,rad)
#' int.circ.line(P1,P1+c(0,1),cent,rad)
#' int.circ.line(P1+c(0,1),P1,cent,rad)
#'
#' dist.point2line(cent,P1,P2)
#' rad2<-dist.point2line(cent,P1,P2)$d
#' int.circ.line(P1,P2,cent,rad2)
#' int.circ.line(P1,P2,cent,rad=.8)
#' int.circ.line(P1,P2,cent,rad=.78)
#'
#' #plot of the line and the circle
#' A<-c(.3,.2); B<-c(.6,.3); cent<-c(1,1); rad<-2 #check dist.point2line(cent,A,B)$dis, .3
#'
#' IPs<-int.circ.line(A,B,cent,rad)
#'
#' xr<-range(A[1],B[1],cent[1])
#' xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
#' x<-seq(xr[1]-rad-xf,xr[2]+rad+xf,l=20)  #try also l=100
#' lnAB<-Line(A,B,x)
#' y<-lnAB$y
#'
#' Xlim<-range(x,cent[1])
#' Ylim<-range(y,A[2],B[2],cent[2]-rad,cent[2]+rad)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(rbind(A,B,cent),pch=1,asp=1,xlab="x",ylab="y",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' lines(x,y,lty=1)
#' interp::circles(cent[1],cent[2],rad)
#' IP.txt<-c()
#' if (!is.null(IPs))
#' {
#'   for (i in 1:(length(IPs)/2))
#'     IP.txt<-c(IP.txt,paste("I",i, sep = ""))
#' }
#' txt<-rbind(A,B,cent,IPs)
#' text(txt+cbind(rep(xd*.03,nrow(txt)),rep(-yd*.03,nrow(txt))),c("A","B","M",IP.txt))
#' }
#'
#' @export int.circ.line
int.circ.line <- function(p1,p2,cent,rad)
{
  if (!is.point(p1) || !is.point(p2) || !is.point(cent) )
  {stop('p1 and p2 and cent must be all numeric 2D points')}

  if (!is.point(rad,1) || rad<0)
  {stop('rad must be a positive scalar')}

  a<-cent[1]; b<-cent[2];r<-rad; c<-p1[1]

  if (p1[1]!=p2[1])
  {
    ln<-Line(p1,p2,0)
    m<-as.numeric(ln$s); d<-as.numeric(ln$int)
    delta<-round(r^2*(1+m^2)-(b-m*a-d)^2,7)
    if (delta<0)
    {int.pts<-NULL
    } else if (delta ==0)
    {
      x<-(a+b*m-d*m)/(1+m^2); y<-(d+a*m+b*m^2)/(1+m^2);
      int.pts<-c(x,y)
    } else
    {
      x1<-(a+b*m-d*m+sqrt(delta))/(1+m^2);
      x2<-(a+b*m-d*m-sqrt(delta))/(1+m^2);
      y1<-(d+a*m+b*m^2+m*sqrt(delta))/(1+m^2);
      y2<-(d+a*m+b*m^2-m*sqrt(delta))/(1+m^2);
      int.pts<-rbind(c(x1,y1),c(x2,y2))

      ord1<-order(int.pts[,1])

      if (p1[1]<p2[1])
      {int.pts<-int.pts[ord1,]
      } else
      {int.pts<-int.pts[rev(ord1),]
      }
    }
  } else
  {
    delta<-round(r^2-(c-a)^2,7)
    if (delta<0)
    {int.pts<-NULL
    } else if (delta ==0)
    {
      x<-c; y<-b;
      int.pts<-c(x,y)
    } else
    {
      x1<-x2<-c;
      y1<-b-sqrt(delta);
      y2<-b+sqrt(delta);
      int.pts<-rbind(c(x1,y1),c(x2,y2))

      ord2<-order(int.pts[,2])

      if (p1[2]<p2[2])
      {int.pts<-int.pts[ord2,]
      } else
      {int.pts<-int.pts[rev(ord2),]
      }

    }
  }
  int.pts
} #end of the function
#'

#################################################################

#' @title The vertices of the Arc Slice (AS) Proximity Region in the standard basic triangle
#'
#' @description Returns the end points of the line segments and arc-slices that constitute the boundary of
#' AS proximity region for a point in the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Vertex regions are based on the center \code{M="CC"} for circumcenter of \eqn{T_b}; or \eqn{M=(m_1,m_2)} in Cartesian
#' coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_b};
#' default is \code{M="CC"} the circumcenter of \eqn{T_b}.
#' \code{rv} is the index of the vertex region \code{p} resides, with default=\code{NULL}.
#'
#' If \code{p} is outside \eqn{T_b}, it returns \code{NULL} for the proximity region.
#' \code{dec} is the number of decimals (default is 4) to round the barycentric coordinates when checking whether
#' the end points fall on the boundary of the triangle \eqn{T_b} or not (so as not to miss the intersection points
#' due to precision in the decimals).
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point whose AS proximity region is to be computed.
#' @param c1,c2 Positive real numbers representing the top vertex in standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \eqn{T_b} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \eqn{T_b}.
#' @param rv The index of the \code{M}-vertex region containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#' @param dec a positive integer the number of decimals (default is 4) to round the barycentric coordinates when checking whether
#' the end points fall on the boundary of the triangle \eqn{T_b} or not.
#'
#' @return A \code{list} with the elements
#' \item{L,R}{The end points of the line segments on the boundary of the AS proximity region.
#'  Each row in \code{L} and \code{R} constitute a line segment on the boundary.}
#' \item{Arc.Slices}{The end points of the arc-slices on the circular parts of the AS proximity region.
#' Here points in row 1 and row 2 constitute the end points of one arc-slice, points on row 3 and row 4
#' constitute the end points for the next arc-slice and so on.}
#'
#' @seealso \code{\link{NAStri}} and \code{\link{IndNASbas.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6  #try also c1<-.2; c2<-.2;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#'
#' set.seed(1)
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' P1<-as.numeric(runif.bas.tri(1,c1,c2)$g);  #try also P1<-c(.3,.2)
#' NASbas.tri(P1,c1,c2)  #default with M="CC"
#' NASbas.tri(P1,c1,c2,M)
#'
#' #or try
#' Rv<-rv.bas.triCC(P1,c1,c2)$rv
#' NASbas.tri(P1,c1,c2,M,Rv)
#'
#' NASbas.tri(c(3,5),c1,c2,M)
#'
#' P2<-c(.5,.4)
#' NASbas.tri(P2,c1,c2,M)
#'
#' P3<-c(1.5,.4)
#' NASbas.tri(P3,c1,c2,M)
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' #plot of the NAS region
#' P1<-as.numeric(runif.bas.tri(1,c1,c2)$g);
#' CC<-circ.cent.bas.tri(c1,c2)
#'
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' rv<-rv.bas.triCC(P1,c1,c2)$rv
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.bas.tri(c1,c2,M)
#' rv<-rv.bas.tri.cent(P1,c1,c2,M)$rv
#' }
#' RV<-Tb[rv,]
#' rad<-Dist(P1,RV)
#'
#' Int.Pts<-NASbas.tri(P1,c1,c2,M)
#'
#' Xlim<-range(Tb[,1],P1[1]+rad,P1[1]-rad)
#' Ylim<-range(Tb[,2],P1[2]+rad,P1[2]-rad)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' points(rbind(Tb,P1,rbind(Int.Pts$L,Int.Pts$R)))
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' interp::circles(P1[1],P1[2],rad,lty=2)
#' L<-Int.Pts$L; R<-Int.Pts$R
#' segments(L[,1], L[,2], R[,1], R[,2], lty=1,col=2)
#' Arcs<-Int.Pts$a;
#' if (!is.null(Arcs))
#' {
#'   K<-nrow(Arcs)/2
#'   for (i in 1:K)
#'   {A1<-Arcs[2*i-1,]; A2<-Arcs[2*i,];
#'   angles<-angle.str2end(A1,P1,A2)$c
#'
#'   plotrix::draw.arc(P1[1],P1[2],rad,angle1=angles[1],angle2=angles[2],col=2)
#'   }
#' }
#'
#' #proximity region with the triangle (i.e., for labeling the vertices of the NAS)
#' IP.txt<-intpts<-c()
#' if (!is.null(Int.Pts$a))
#' {
#'  intpts<-unique(round(Int.Pts$a,7))
#'   #this part is for labeling the intersection points of the spherical
#'   for (i in 1:(length(intpts)/2))
#'     IP.txt<-c(IP.txt,paste("I",i+1, sep = ""))
#' }
#' txt<-rbind(Tb,P1,cent,intpts)
#' txt.str<-c("A","B","C","P1",cent.name,IP.txt)
#' text(txt+cbind(rep(xd*.02,nrow(txt)),rep(-xd*.03,nrow(txt))),txt.str)
#'
#' c1<-.4; c2<-.6;
#' P1<-c(.3,.2)
#' NASbas.tri(P1,c1,c2,M)
#' }
#'
#' @export NASbas.tri
NASbas.tri <- function(p,c1,c2,M="CC",rv=NULL,dec=4)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

    if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  CC = circ.cent.tri(Tb)
  if (identical(M,"CC") )
    { M<-CC }

    if (dimension(M)==3)
    {M<-bary2cart(M,Tb)}

    if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,Tb,boundary=FALSE)$in.tri))
    {stop('center is not the circumcenter or not in the interior of the triangle')}

  #If p is outside the closed triangle
  if (!in.triangle(p,Tb,boundary=TRUE)$in.tri)
  {reg<-list(L=NULL, R=NULL, Arc.Slices=NULL); return(reg); stop}

  #If p is at a vertex of the triangle
  if (dist.point2set(p,Tb)$dist==0)
  {reg<-list(L=p, R=p, Arc.Slices=NULL); return(reg); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.bas.triCC(p,c1,c2)$rv,rv.bas.tri.cent(p,c1,c2,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}


  Seg.LPts<-Seg.RPts<-Arc.Pts<-NULL; #segment and arc end points
  if (rv==1)
  {
    rad<-Dist(p,y1)
    if (in.circle(y1,p,rad,boundary = TRUE) && in.circle(y2,p,rad,boundary = TRUE) && in.circle(y3,p,rad,boundary = TRUE))
    {
      Seg.LPts<-Tb; Seg.RPts<-rbind(y2,y3,y1)
      Arc.Pts<-NULL
    } else
    {
      pts1<-int.circ.line(y1,y2,p,rad); pts2<-int.circ.line(y2,y3,p,rad); pts3<-int.circ.line(y1,y3,p,rad)

      if (length(pts1)/2>1)
      {p1<-y1; p2<-pts1[2,]
      cond<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0);
      if (cond)
      {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
      Arc.Pts<-rbind(Arc.Pts,p2)
      } else
      {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y2)
      }
      } else
      {Arc.Pts<-rbind(Arc.Pts,pts1,pts1)}

      if (length(pts2)/2>1)
      {p1<-pts2[1,]; p2<-pts2[2,];
      cond1<-all(round(in.triangle(p1,Tb,boundary = TRUE)$b,dec)>=0);
      cond2<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0)
      if (cond1 && cond2)
      {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
      Arc.Pts<-rbind(Arc.Pts,p1,p2)
      } else
      { if (cond1) {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y3); Arc.Pts<-rbind(Arc.Pts,p1)}
        if (cond2) {Seg.LPts<-rbind(Seg.LPts,p2); Seg.RPts<-rbind(Seg.RPts,y2);Arc.Pts<-rbind(Arc.Pts,p2)}
      }
      } else
      {Arc.Pts<-rbind(Arc.Pts,pts2,pts2)}

      if (length(pts3)/2>1)
      {p1<-y1; p2<-pts3[2,]
      cond<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0);
      if (cond)
      {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
      Arc.Pts<-rbind(Arc.Pts,p2)
      } else
      {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y3)
      }
      } else
      {Arc.Pts<-rbind(Arc.Pts,pts3,pts3)}
    }
  } else
    if (rv==2)
    {
      rad<-Dist(p,y2)
      if (in.circle(y1,p,rad,boundary = TRUE) && in.circle(y2,p,rad,boundary = TRUE) && in.circle(y3,p,rad,boundary = TRUE))
      {
        Seg.LPts<-Tb; Seg.RPts<-rbind(y2,y3,y1)
        Arc.Pts<-NULL
      } else
      {
        pts1<-int.circ.line(y2,y3,p,rad); pts2<-int.circ.line(y3,y1,p,rad); pts3<-int.circ.line(y2,y1,p,rad)

        if (length(pts1)/2>1)
        {
          p1<-y2; p2<-pts1[2,]
          cond<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0);
          if (cond)
          {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
          Arc.Pts<-rbind(Arc.Pts,p2)
          } else
          {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y3)
          }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts1,pts1)}

        if (length(pts2)/2>1)
        { p1<-pts2[1,]; p2<-pts2[2,];
        cond1<-all(round(in.triangle(p1,Tb,boundary = TRUE)$b,dec)>=0);
        cond2<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0)
        if (cond1 && cond2)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p1,p2)
        } else
        {if (cond1) {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y1); Arc.Pts<-rbind(Arc.Pts,p1)}
          if (cond2) {Seg.LPts<-rbind(Seg.LPts,p2); Seg.RPts<-rbind(Seg.RPts,y3); Arc.Pts<-rbind(Arc.Pts,p2)}
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts2,pts2)}

        if (length(pts3)/2>1)
        {p1<-y2; p2<-pts3[2,]
        cond<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0);
        if (cond)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p2)
        } else
        { Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y1)
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts3,pts3)}
      }
    } else {
      rad<-Dist(p,y3)
      if (in.circle(y1,p,rad,boundary = TRUE) && in.circle(y2,p,rad,boundary = TRUE) && in.circle(y3,p,rad,boundary = TRUE))
      {
        Seg.LPts<-Tb; Seg.RPts<-rbind(y2,y3,y1)
        Arc.Pts<-NULL
      } else
      {
        pts1<-int.circ.line(y3,y1,p,rad); pts2<-int.circ.line(y1,y2,p,rad); pts3<-int.circ.line(y3,y2,p,rad)

        if (length(pts1)/2>1)
        { p1<-y3; p2<-pts1[2,]
        cond<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0);
        if (cond)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p2)
        } else
        {
          Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y1)
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts1,pts1)}

        if (length(pts2)/2>1)
        {p1<-pts2[1,]; p2<-pts2[2,];
        cond1<-all(round(in.triangle(p1,Tb,boundary = TRUE)$b,dec)>=0);
        cond2<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0)
        if (cond1 && cond2)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p1,p2)
        } else
        {if (cond1) {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y2); Arc.Pts<-rbind(Arc.Pts,p1)}
          if (cond2) {Seg.LPts<-rbind(Seg.LPts,y1); Seg.RPts<-rbind(Seg.RPts,p2); Arc.Pts<-rbind(Arc.Pts,p2)}
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts2,pts2)}

        if (length(pts3)/2>1)
        {p1<-y3; p2<-pts3[2,]
        cond<-all(round(in.triangle(p2,Tb,boundary = TRUE)$b,dec)>=0);
        if (cond)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p2)
        } else
        {
          Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y2)
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts3,pts3)}
      }
    }

  row.names(Seg.LPts) <- row.names(Seg.RPts) <- row.names(Arc.Pts) <- NULL

  #to determine the angles between the vectors joining arc slice end points to the point p
  #and the horizontal line crossing p
  Angles = NULL
  if (!is.null(Arc.Pts))
  {
    K<-nrow(Arc.Pts)
    for (i in 1:K)
    { vec<-Arc.Pts[i,]-p
    Angles<-c(Angles,acos(vec[1]/sqrt(sum(vec^2))) )
    }
  }

  list(L=Seg.LPts,
       R=Seg.RPts,
       arc.slices=Arc.Pts,
       Angles=Angles)
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Arc Slice Proximity Catch Digraphs
#' (AS-PCDs) - standard basic triangle case
#'
#' @description Returns \eqn{I(p2 \in N_{AS}(p1))} for points \code{p1} and \code{p2}, that is, returns 1 if \eqn{p2} is in \eqn{N_{AS}(p1)}, returns 0
#' otherwise, where \eqn{N_{AS}(x)} is the AS proximity region for point \eqn{x}.
#'
#' AS proximity region is constructed in the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Vertex regions are based on the center \code{M="CC"} for circumcenter of \eqn{T_b};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of \eqn{T_b}; default is \code{M="CC"} i.e., circumcenter of \eqn{T_b}.
#' \code{rv} is the index of the vertex region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_b}, the function returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p1 A 2D point whose AS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the AS proximity region of
#' \code{p1} or not.
#' @param c1,c2 Positive real numbers representing the top vertex in standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \eqn{T_b}.
#' @param rv The index of the \code{M}-vertex region in \eqn{T_b} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(p2 \in N_{AS}(p1))} for points \code{p1} and \code{p2}, that is, returns 1 if \eqn{p2} is in \eqn{N_{AS}(p1)}
#' (i.e., if there is an arc from \code{p1} to \code{p2}), returns 0 otherwise.
#'
#' @seealso \code{\link{IndNAStri}} and \code{\link{NAStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' P1<-as.numeric(runif.bas.tri(1,c1,c2)$g)
#' P2<-as.numeric(runif.bas.tri(1,c1,c2)$g)
#' IndNASbas.tri(P1,P2,c1,c2,M)
#'
#' P1<-c(.3,.2)
#' P2<-c(.6,.2)
#' IndNASbas.tri(P1,P2,c1,c2,M)
#'
#' #or try
#' Rv<-rv.bas.triCC(P1,c1,c2)$rv
#' IndNASbas.tri(P1,P2,c1,c2,M,Rv)
#'
#' P1<-c(.3,.2)
#' P2<-c(.8,.2)
#' IndNASbas.tri(P1,P2,c1,c2,M)
#'
#' P3<-c(.5,.4)
#' IndNASbas.tri(P1,P3,c1,c2,M)
#'
#' P4<-c(1.5,.4)
#' IndNASbas.tri(P1,P4,c1,c2,M)
#' IndNASbas.tri(P4,P4,c1,c2,M)
#'
#' c1<-.4; c2<-.6;
#' P1<-c(.3,.2)
#' P2<-c(.6,.2)
#' IndNASbas.tri(P1,P2,c1,c2,M)
#' }
#'
#' @export IndNASbas.tri
IndNASbas.tri <- function(p1,p2,c1,c2,M="CC",rv=NULL)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  CC = circ.cent.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,Tb,boundary=TRUE)$in.tri || !in.triangle(p2,Tb,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.bas.triCC(p1,c1,c2)$rv,rv.bas.tri.cent(p1,c1,c2,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  X1<-p1[1]; Y1<-p1[2];
  X2<-p2[1]; Y2<-p2[2];
  arc<-0;
  if (rv==1)
  {
    if ( (Y2-Y1)^2+(X2-X1)^2 <= X1^2+Y1^2 ) {arc <-1}
  } else {
    if (rv==2)
    {
      if ( (Y2-Y1)^2+(X2-X1)^2 <= (X1-1)^2+Y1^2 ) {arc <-1}
    } else {
      if ( (Y2-Y1)^2+(X2-X1)^2 <= (X1-c1)^2+(Y1-c2)^2 ) {arc<-1}
    }}
  arc
} #end of the function
#'

#################################################################

#' @title The vertices of the Arc Slice (AS) Proximity Region in a general triangle
#'
#' @description Returns the end points of the line segments and arc-slices that constitute the boundary of AS proximity
#' region for a point in the triangle \code{tri}\eqn{=T(A,B,C)=}\code{(rv=1,rv=2,rv=3)}.
#'
#' Vertex regions are based on the center \code{M="CC"} for circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \code{M="CC"} the circumcenter of \code{tri}. \code{rv} is the index of the vertex region \code{p1} resides,
#' with default=\code{NULL}.
#'
#' If \code{p} is outside of \code{tri}, it returns \code{NULL} for the proximity region.
#' \code{dec} is the number of decimals (default is 4) to round the barycentric coordinates when checking the points
#' fall on the boundary of the triangle \code{tri} or not (so as not to miss the intersection points due to precision
#' in the decimals).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point whose AS proximity region is to be computed.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates which serves as a center in the interior of the triangle \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param rv Index of the \code{M}-vertex region containing the point \code{p}, either \code{1,2,3} or
#' \code{NULL} (default is \code{NULL}).
#' @param dec a positive integer the number of decimals (default is 4) to round the
#' barycentric coordinates when checking whether
#' the end points fall on the boundary of the triangle \code{tri} or not.
#'
#' @return A \code{list} with the elements
#' \item{L,R}{End points of the line segments on the boundary of the AS proximity region.
#'  Each row in \code{L} and \code{R} constitute a pair of points that determine a line segment on the boundary.}
#' \item{arc.slices}{The end points of the arc-slices on the circular parts of the AS proximity region.
#' Here points in rows 1 and 2 constitute the end points of the first arc-slice, points on rows 3 and 4
#' constitute the end points for the next arc-slice and so on.}
#' \item{Angles}{The angles (in radians) between the vectors joining arc slice end points to the point \code{p}
#' with the horizontal line crossing the point \code{p}}
#'
#' @seealso \code{\link{NASbas.tri}}, \code{\link{NPEtri}}, \code{\link{NCStri}} and \code{\link{IndNAStri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(.6,.2)
#'
#' P1<-as.numeric(runif.tri(1,Tr)$g)  #try also P1<-c(1.3,1.2)
#' NAStri(P1,Tr,M)
#'
#' #or try
#' Rv<-rv.triCC(P1,Tr)$rv
#' NAStri(P1,Tr,M,Rv)
#'
#' NAStri(c(3,5),Tr,M)
#'
#' P2<-c(1.5,1.4)
#' NAStri(P2,Tr,M)
#'
#' P3<-c(1.5,.4)
#' NAStri(P3,Tr,M)
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#'
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' rv<-rv.triCC(P1,Tr)$rv
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.tri(Tr,M)
#' rv<-rv.tri.cent(P1,Tr,M)$rv
#' }
#' RV<-Tr[rv,]
#' rad<-Dist(P1,RV)
#'
#' Int.Pts<-NAStri(P1,Tr,M)
#'
#' #plot of the NAS region
#' Xlim<-range(Tr[,1],P1[1]+rad,P1[1]-rad)
#' Ylim<-range(Tr[,2],P1[2]+rad,P1[2]-rad)
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' #asp=1 must be the case to have the arc properly placed in the figure
#' polygon(Tr)
#' points(rbind(Tr,P1,rbind(Int.Pts$L,Int.Pts$R)))
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' interp::circles(P1[1],P1[2],rad,lty=2)
#' L<-Int.Pts$L; R<-Int.Pts$R
#' segments(L[,1], L[,2], R[,1], R[,2], lty=1,col=2)
#' Arcs<-Int.Pts$a;
#' if (!is.null(Arcs))
#' {
#'   K<-nrow(Arcs)/2
#'   for (i in 1:K)
#'   {A1<-Int.Pts$arc[2*i-1,]; A2<-Int.Pts$arc[2*i,];
#'   angles<-angle.str2end(A1,P1,A2)$c
#'
#'   test.ang1<-angles[1]+(.01)*(angles[2]-angles[1])
#'   test.Pnt<-P1+rad*c(cos(test.ang1),sin(test.ang1))
#'   if (!in.triangle(test.Pnt,Tr,boundary = TRUE)$i) {angles<-c(min(angles),max(angles)-2*pi)}
#'   plotrix::draw.arc(P1[1],P1[2],rad,angle1=angles[1],angle2=angles[2],col=2)
#'   }
#' }
#'
#' #proximity region with the triangle (i.e., for labeling the vertices of the NAS)
#' IP.txt<-intpts<-c()
#' if (!is.null(Int.Pts$a))
#' {
#'  intpts<-unique(round(Int.Pts$a,7))
#'   #this part is for labeling the intersection points of the spherical
#'   for (i in 1:(length(intpts)/2))
#'     IP.txt<-c(IP.txt,paste("I",i+1, sep = ""))
#' }
#' txt<-rbind(Tr,P1,cent,intpts)
#' txt.str<-c("A","B","C","P1",cent.name,IP.txt)
#' text(txt+cbind(rep(xd*.02,nrow(txt)),rep(-xd*.03,nrow(txt))),txt.str)
#'
#' P1<-c(.3,.2)
#' NAStri(P1,Tr,M)
#' }
#'
#' @export NAStri
NAStri <- function(p,tri,M="CC",rv=NULL,dec=4)
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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  #If p is outside the closed triangle
  if (!in.triangle(p,tri,boundary=TRUE)$in.tri)
  {reg<-list(L=NULL, R=NULL, Arc.Slices=NULL); return(reg); stop}

  #If p is at a vertex of the triangle
  if (dist.point2set(p,tri)$dist==0)
  {reg<-list(L=p, R=p, Arc.Slices=NULL); return(reg); stop}

  Tr<-tri[order(tri[,1]),] #order the vertices according to their x axis, so that angles for the
  ifelse(Tr[2,2]>Tr[3,2],tri<-Tr[c(1,3,2),],tri<-Tr)   #arcs are provided in counter-clockwise

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p,tri)$rv,rv.tri.cent(p,tri,M)$rv)  #vertex region for p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];

  Seg.LPts<-Seg.RPts<-Arc.Pts<-NULL; #segment and arc end points
  if (rv==1)
  { rad<-Dist(p,y1)
  if (in.circle(y1,p,rad,boundary = TRUE) && in.circle(y2,p,rad,boundary = TRUE) && in.circle(y3,p,rad,boundary = TRUE))
  {
    Seg.LPts<-tri; Seg.RPts<-rbind(y2,y3,y1)
    Arc.Pts<-NULL
  } else
  {
    pts1<-int.circ.line(y1,y2,p,rad); pts2<-int.circ.line(y2,y3,p,rad); pts3<-int.circ.line(y1,y3,p,rad)

    if (length(pts1)/2>1)
    {p1<-y1; p2<-pts1[2,]
    cond<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0);
    if (cond)
    {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
    Arc.Pts<-rbind(Arc.Pts,p2)
    } else
    {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y2)
    }
    } else
    {Arc.Pts<-rbind(Arc.Pts,pts1,pts1)}

    if (length(pts2)/2>1)
    {p1<-pts2[1,]; p2<-pts2[2,];
    cond1<-all(round(in.triangle(p1,tri,boundary = TRUE)$b,dec)>=0);
    cond2<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0)
    if (cond1 && cond2)
    {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
    Arc.Pts<-rbind(Arc.Pts,p1,p2)
    } else
    { if (cond1) {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y3); Arc.Pts<-rbind(Arc.Pts,p1)}
      if (cond2) {Seg.LPts<-rbind(Seg.LPts,p2); Seg.RPts<-rbind(Seg.RPts,y2); Arc.Pts<-rbind(Arc.Pts,p2)}
    }
    } else
    {Arc.Pts<-rbind(Arc.Pts,pts2,pts2)}

    if (length(pts3)/2>1)
    {p1<-y1; p2<-pts3[2,]
    cond<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0);
    if (cond)
    {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
    Arc.Pts<-rbind(Arc.Pts,p2)
    } else
    {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y3)
    }
    } else
    {Arc.Pts<-rbind(Arc.Pts,pts3,pts3)}
  }
  } else
    if (rv==2)
    {
      rad<-Dist(p,y2)
      if (in.circle(y1,p,rad,boundary = TRUE) && in.circle(y2,p,rad,boundary = TRUE) && in.circle(y3,p,rad,boundary = TRUE))
      {
        Seg.LPts<-tri; Seg.RPts<-rbind(y2,y3,y1)
        Arc.Pts<-NULL
      } else
      {
        pts1<-int.circ.line(y2,y3,p,rad); pts2<-int.circ.line(y3,y1,p,rad); pts3<-int.circ.line(y2,y1,p,rad)

        if (length(pts1)/2>1)
        {
          p1<-y2; p2<-pts1[2,]
          cond<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0);
          if (cond)
          {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
          Arc.Pts<-rbind(Arc.Pts,p2)
          } else
          {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y3)
          }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts1,pts1)}

        if (length(pts2)/2>1)
        { p1<-pts2[1,]; p2<-pts2[2,];
        cond1<-all(round(in.triangle(p1,tri,boundary = TRUE)$b,dec)>=0);
        cond2<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0)
        if (cond1 && cond2)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p1,p2)
        } else
        {if (cond1) {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y1); Arc.Pts<-rbind(Arc.Pts,p1)}
          if (cond2) {Seg.LPts<-rbind(Seg.LPts,p2); Seg.RPts<-rbind(Seg.RPts,y3); Arc.Pts<-rbind(Arc.Pts,p2)}
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts2,pts2)}

        if (length(pts3)/2>1)
        {p1<-y2; p2<-pts3[2,]
        cond<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0);
        if (cond)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p2)
        } else
        { Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y1)
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts3,pts3)}
      }
    } else {
      rad<-Dist(p,y3)
      if (in.circle(y1,p,rad,boundary = TRUE) && in.circle(y2,p,rad,boundary = TRUE) && in.circle(y3,p,rad,boundary = TRUE))
      {
        Seg.LPts<-tri; Seg.RPts<-rbind(y2,y3,y1)
        Arc.Pts<-NULL
      } else
      {
        pts1<-int.circ.line(y3,y1,p,rad); pts2<-int.circ.line(y1,y2,p,rad); pts3<-int.circ.line(y3,y2,p,rad)

        if (length(pts1)/2>1)
        { p1<-y3; p2<-pts1[2,]
        cond<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0);
        if (cond)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p2)
        } else
        {
          Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y1)
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts1,pts1)}

        if (length(pts2)/2>1)
        {p1<-pts2[1,]; p2<-pts2[2,];
        cond1<-all(round(in.triangle(p1,tri,boundary = TRUE)$b,dec)>=0);
        cond2<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0)
        if (cond1 && cond2)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p1,p2)
        } else
        {if (cond1) {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y2); Arc.Pts<-rbind(Arc.Pts,p1)}
          if (cond2) {Seg.LPts<-rbind(Seg.LPts,y1); Seg.RPts<-rbind(Seg.RPts,p2); Arc.Pts<-rbind(Arc.Pts,p2)}
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts2,pts2)}

        if (length(pts3)/2>1)
        {p1<-y3; p2<-pts3[2,]
        cond<-all(round(in.triangle(p2,tri,boundary = TRUE)$b,dec)>=0);
        if (cond)
        {Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,p2)
        Arc.Pts<-rbind(Arc.Pts,p2)
        } else
        {
          Seg.LPts<-rbind(Seg.LPts,p1); Seg.RPts<-rbind(Seg.RPts,y2)
        }
        } else
        {Arc.Pts<-rbind(Arc.Pts,pts3,pts3)}
      }
    }

  rownames(Seg.LPts)<-rownames(Seg.RPts)<-rownames(Arc.Pts)<-NULL #removing the row names

  #to determine the angles between the vectors joining arc slice end points to the point p
  #and the horizontal line crossing p
  Angles = NULL
  if (!is.null(Arc.Pts))
  {
    K<-nrow(Arc.Pts)
    for (i in 1:K)
    { vec<-Arc.Pts[i,]-p
    Angles<-c(Angles,acos(vec[1]/sqrt(sum(vec^2))) )
    }
  }

  list(L=Seg.LPts,
       R=Seg.RPts,
       arc.slices=Arc.Pts,
       Angles=Angles)
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Arc Slice Proximity Catch Digraphs
#' (AS-PCDs) - one triangle case
#'
#' @description Returns \eqn{I(p2 \in N_{AS}(p1))} for points \code{p1} and \code{p2}, that is, returns 1 if \eqn{p2} is in \eqn{N_{AS}(p1)}, returns 0
#' otherwise, where \eqn{N_{AS}(x)} is the AS proximity region for point \eqn{x}.
#'
#' AS proximity regions are constructed with respect to the triangle, \code{tri}\eqn{=T(A,B,C)=}\code{(rv=1,rv=2,rv=3)},
#' and vertex regions are based on the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
#' \code{rv} is the index of the vertex region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \code{tri}, the function returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p1 A 2D point whose AS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the AS proximity region of
#' \code{p1} or not.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param rv The index of the \code{M}-vertex region in \code{tri} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(p2 \in N_{AS}(p1))} for \code{p1}, that is, returns 1 if \eqn{p2} is in \eqn{N_{AS}(p1)}, returns 0 otherwise
#'
#' @seealso \code{\link{IndNASbas.tri}}, \code{\link{IndNPEtri}}, and \code{\link{IndNCStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#'
#' Tr<-rbind(A,B,C);
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' P1<-as.numeric(runif.tri(1,Tr)$g)
#' P2<-as.numeric(runif.tri(1,Tr)$g)
#' IndNAStri(P1,P2,Tr,M)
#'
#' P1<-c(1.3,1.2)
#' P2<-c(1.8,.5)
#' IndNAStri(P1,P2,Tr,M)
#' IndNAStri(P1,P1,Tr,M)
#'
#' #or try
#' Rv<-rv.triCC(P1,Tr)$rv
#' IndNAStri(P1,P2,Tr,M,Rv)
#'
#' P3<-c(1.6,1.4)
#' IndNAStri(P1,P3,Tr,M)
#'
#' P4<-c(1.5,1.0)
#' IndNAStri(P1,P4,Tr,M)
#'
#' P5<-c(.5,1.0)
#' IndNAStri(P1,P5,Tr,M)
#' IndNAStri(P5,P5,Tr,M)
#'
#' #or try
#' Rv<-rv.triCC(P5,Tr)$rv
#' IndNAStri(P5,P5,Tr,M,Rv)
#' }
#'
#' @export IndNAStri
IndNAStri <- function(p1,p2,tri,M="CC",rv=NULL)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,tri,boundary=TRUE)$in.tri || !in.triangle(p2,tri,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p1,tri)$rv,rv.tri.cent(p1,tri,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  X1<-p1[1]; Y1<-p1[2]; X2<-p2[1]; Y2<-p2[2];
  a1<-tri[1,1]; a2<-tri[1,2]; b1<-tri[2,1]; b2<-tri[2,2]; c1<-tri[3,1]; c2<-tri[3,2];
  arc<-0;
  if (rv==1)
  {
    if ( (Y2-Y1)^2+(X2-X1)^2 <= (X1-a1)^2+(Y1-a2)^2 ) {arc <-1}
  } else {
    if (rv==2)
    {
      if ( (Y2-Y1)^2+(X2-X1)^2 <= (X1-b1)^2+(Y1-b2)^2 ) {arc <-1}
    } else {
      if ( (Y2-Y1)^2+(X2-X1)^2 <= (X1-c1)^2+(Y1-c2)^2 ) {arc<-1}
    }}
  arc
} #end of the function
#'

#################################################################

#' @title Number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case
#'
#' @description Returns the number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs)
#' whose vertices are the 2D data set, \code{Xp}. It also provides number of vertices (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' The data points could be inside or outside a general
#' triangle \code{tri}\eqn{=T(A,B,C)=}\code{(rv=1,rv=2,rv=3)}, with vertices of \code{tri} stacked row-wise.
#'
#' AS proximity regions are defined with respect to the triangle \code{tri} and vertex regions are
#' based on the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
#' For the number of arcs, loops are not allowed,
#' so arcs are only possible for points inside the triangle, \code{tri}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the digraph (i.e., AS-PCD).
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{num.arcs}{Number of arcs of the AS-PCD}
#' \item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points that reside in the triangle}
#'
#' @seealso \code{\link{NumArcsAS}}, \code{\link{NumArcsPEtri}}, and \code{\link{NumArcsCStri}}
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
#' n<-10  #try also n<-20
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' NumArcsAStri(Xp,Tr,M)
#'
#' Xp2=Xp+1
#' NumArcsAStri(Xp2,Tr,M)
#'
#' NumArcsAStri(rbind(Xp,c(0,2)),Tr,M)
#' }
#'
#' @export NumArcsAStri
NumArcsAStri <- function(Xp,tri,M="CC")
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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  n<-nrow(Xp)
  arcs<-0
  ind.in.tri = NULL
  if (n<=0)
  {
    arcs<-0
  } else
  {
    for (i in 1:n)
    { p1<-Xp[i,]
    if (in.triangle(p1,tri,boundary=TRUE)$in.tri)
    { vert<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p1,tri)$rv,rv.tri.cent(p1,tri,M)$rv)  #vertex region for pt
    ind.in.tri = c(ind.in.tri,i)

     for (j in (1:n)[-i])  #to avoid loops
    {
      arcs<-arcs+IndNAStri(p1,Xp[j,],tri,M,rv=vert)
    }
    }
    }
  }

  NinTri = length(ind.in.tri)

  res<-list(num.arcs=arcs, #number of arcs for the AS-PCD
            num.in.tri=NinTri, # number of Xp points in CH of Yp points
            ind.in.tri=ind.in.tri) #indices of data points inside the triangle

res
} #end of the function
#'

#################################################################

#' @title Arc density of Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case
#'
#' @description Returns the arc density of AS-PCD whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the triangle \code{tri}.
#'
#' AS proximity regions is defined with respect to \code{tri}
#' and vertex regions are defined with the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
#' For the number of arcs, loops are not allowed so arcs are only possible for points inside \code{tri}
#' for this function.
#'
#' \code{tri.cor} is a logical argument for triangle correction (default is \code{TRUE}),
#' if \code{TRUE}, only the points inside the triangle are considered (i.e., digraph induced by these vertices
#' are considered) in computing the arc density, otherwise all points are considered
#' (for the number of vertices in the denominator of arc density).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param tri.cor A logical argument for computing the arc density for only the points inside the triangle, \code{tri}
#' (default is \code{tri.cor=FALSE}), i.e., if \code{tri.cor=TRUE} only the induced digraph with the vertices inside \code{tri} are considered in the
#' computation of arc density.
#'
#' @return Arc density of AS-PCD whose vertices are the 2D numerical data set, \code{Xp};
#' AS proximity regions are defined with respect to the triangle \code{tri} and \eqn{CC}-vertex regions.
#'
#' @seealso \code{\link{ASarcdens.tri}}, \code{\link{CSarcdens.tri}}, and \code{\link{NumArcsAStri}}
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
#' set.seed(1)
#' n<-10  #try also n<-20
#'
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' NumArcsAStri(Xp,Tr,M)
#' ASarcdens.tri(Xp,Tr,M)
#' ASarcdens.tri(Xp,Tr,M,tri.cor = FALSE)
#'
#' ASarcdens.tri(Xp,Tr,M)
#' }
#'
#' @export ASarcdens.tri
ASarcdens.tri <- function(Xp,tri,M="CC",tri.cor=FALSE)
{
  if (!is.numeric(as.matrix(Xp)) )
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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  nx<-nrow(Xp)
  narcs<-NumArcsAStri(Xp,tri,M)$num.arcs

  if (tri.cor==TRUE)
  {
    ind.it<-c()
    for (i in 1:nx)
    {
      ind.it<-c(ind.it,in.triangle(Xp[i,],tri)$in.tri)
    }
    dat.it<-Xp[ind.it,] #Xp points inside the triangle
    NinTri<-nrow(dat.it)
    if (NinTri<=1)
    {stop('not enough points in the triangle to compute the arc density')}
    n<-NinTri
  } else
  {
    n<-nx
  }
  rho<-narcs/(n*(n-1))
  rho
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point in set \code{S} to the point \code{p} for
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case
#'
#' @description Returns I(\eqn{pt \in N_{AS}(x)} for some \eqn{x \in S}), that is, returns 1 if \eqn{p} is in \eqn{\cup_{x \in S}N_{AS}(x)},
#' returns 0 otherwise, where \eqn{N_{AS}(x)} is the AS proximity region for point \eqn{x}.
#'
#' AS proximity regions are constructed with respect to the triangle, \code{tri}\eqn{=T(A,B,C)=}\code{(rv=1,rv=2,rv=3)},
#' and vertices of \code{tri} are also labeled as 1,2, and 3, respectively.
#'
#' Vertex regions are based on the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
#'
#' If \code{p} is not in \code{S} and either \code{p} or all points in \code{S} are outside \code{tri}, it returns 0,
#' but if \code{p} is in \code{S}, then it always returns 1 (i.e., loops are allowed).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points whose AS proximity regions are considered.
#' @param p A 2D point. The function determines whether \code{p} is inside the union of AS proximity
#' regions of points in \code{S} or not.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#'
#' @return \eqn{I(pt \in \cup_{x in S}N_{AS}(x,r))}, that is, returns 1 if \code{p} is in \code{S} or inside \eqn{N_{AS}(x)} for at least
#' one \eqn{x} in \code{S}, returns 0 otherwise, where AS proximity region is constructed in \code{tri}
#'
#' @seealso \code{\link{IndNAStri}}, \code{\link{IndNAStriSet}}, and \code{\link{IndNCStriSet}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$gen.points
#'
#' S<-rbind(Xp[1,],Xp[2,])  #try also S<-c(1.5,1)
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' IndNAStriSet(S,Xp[3,],Tr,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' IndNAStriSet(S,Xp[3,],Tr,M)
#'
#' IndNAStriSet(S,Xp[6,],Tr,M)
#'
#' S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
#' IndNAStriSet(S,Xp[3,],Tr,M)
#'
#' IndNAStriSet(c(.2,.5),Xp[2,],Tr,M)
#' IndNAStriSet(Xp,c(.2,.5),Tr,M)
#' IndNAStriSet(Xp,Xp[2,],Tr,M)
#' IndNAStriSet(c(.2,.5),c(.2,.5),Tr,M)
#' IndNAStriSet(Xp[5,],Xp[2,],Tr,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,],c(.2,.5))
#' IndNAStriSet(S,Xp[3,],Tr,M)
#'
#' P<-c(.4,.2)
#' S<-Xp[c(1,3,4),]
#' IndNAStriSet(Xp,P,Tr,M)
#' IndNAStriSet(S,P,Tr,M)
#'
#' IndNAStriSet(rbind(S,S),P,Tr,M)
#' }
#'
#' @export IndNAStriSet
IndNAStriSet <- function(S,p,tri,M="CC")
{
  if (!is.numeric(as.matrix(S)))
  {stop('S must be a matrix of numeric values')}

  if (is.point(S))
  { S<-matrix(S,ncol=2)
  } else
  {S<-as.matrix(S)
  if (ncol(S)!=2 )
  {stop('S must be of dimension nx2')}
  }

  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  k<-nrow(S);
  dom<-0; i<-1;
  while (dom ==0 && i<= k)
  {
    if (IndNAStri(S[i,],p,tri,M)==1)
    {dom<-1};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the set of points \code{S} being a dominating set or not for Arc Slice Proximity
#' Catch Digraphs (AS-PCDs) - one triangle case
#'
#' @description Returns \eqn{I(}\code{S} a dominating set of AS-PCD\eqn{)}, that is, returns 1 if \code{S} is a dominating set of AS-PCD,
#' returns 0 otherwise.
#'
#' AS-PCD has vertex set \code{Xp} and AS proximity region is constructed with vertex
#' regions based on the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}
#' whose vertices are also labeled as edges 1, 2, and 3, respectively.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points which is to be tested for being a dominating set for the AS-PCDs.
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#'
#' @return \eqn{I(}\code{S} a dominating set of AS-PCD\eqn{)}, that is, returns 1 if \code{S} is a dominating set of AS-PCD whose
#' vertices are the data points in \code{Xp}; returns 0 otherwise, where AS proximity region is constructed in
#' the triangle \code{tri}.
#'
#' @seealso \code{\link{IndNAStriSet}}, \code{\link{IndNPEtri.domset}} and \code{\link{IndNCStri.domset}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#'
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$gen.points
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' S<-rbind(Xp[1,],Xp[2,])
#' IndNAStri.domset(S,Xp,Tr,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' IndNAStri.domset(S,Xp,Tr,M)
#'
#' S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
#' IndNAStri.domset(S,Xp,Tr,M)
#'
#' IndNAStri.domset(c(.2,.5),Xp,Tr,M)
#' IndNAStri.domset(c(.2,.5),c(.2,.5),Tr,M)
#' IndNAStri.domset(Xp[5,],Xp[2,],Tr,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,],c(.2,.5))
#' IndNAStri.domset(S,Xp[3,],Tr,M)
#'
#' IndNAStri.domset(Xp,Xp,Tr,M)
#'
#' P<-c(.4,.2)
#' S<-Xp[c(1,3,4),]
#' IndNAStri.domset(Xp,P,Tr,M)
#' IndNAStri.domset(S,P,Tr,M)
#' IndNAStri.domset(S,Xp,Tr,M)
#'
#' IndNAStri.domset(rbind(S,S),Xp,Tr,M)
#' }
#'
#' @export IndNAStri.domset
IndNAStri.domset <- function(S,Xp,tri,M="CC")
{
  if (!is.numeric(as.matrix(S)) || !is.numeric(as.matrix(Xp)))
  {stop('S and Xp must be numeric')}

  if (is.point(S))
  { S<-matrix(S,ncol=2)
  } else
  {S<-as.matrix(S)
  if (ncol(S)!=2 )
  {stop('S must be of dimension nx2')}
  }

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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  k<-nrow(S);
  n<-nrow(Xp);

  dom<-1; i<-1;
  while (dom ==1 && i<= n)
  {
    if (IndNAStriSet(S,Xp[i,],tri,M)==0)  #this is where tri is used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title Indicator for an upper bound for the domination number of Arc Slice Proximity Catch Digraph
#' (AS-PCD) by the exact algorithm - one triangle case
#'
#' @description Returns \eqn{I(}domination number of AS-PCD whose vertices are the data points \code{Xp} is less than or equal to \code{k}\eqn{)},
#' that is, returns 1 if the domination number of AS-PCD is less than the prespecified value \code{k}, returns 0
#' otherwise. It also provides the vertices (i.e., data points) in a dominating set of size \code{k} of AS-PCD.
#'
#' AS proximity regions are constructed with respect to the triangle \code{tri} and
#' vertex regions are based on the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
#'
#' The vertices of \code{tri}angle, \code{tri},
#' are labeled as \eqn{1,2,3} according to the row number the vertex is recorded in \code{tri}.
#' Loops are allowed in the digraph.
#' It takes a long time for large number of vertices (i.e., large number of row numbers).
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param k A positive integer to be tested for an upper bound for the domination number of AS-PCDs.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{domUB}{The suggested upper bound (to be checked) for the domination number of AS-PCD.
#' It is prespecified as \code{k} in the function arguments.}
#' \item{IndUBdom}{The indicator for the upper bound for domination number of AS-PCD being the
#' specified value \code{k} or not. It returns 1 if the upper bound is \code{k}, and 0 otherwise.}
#' \item{ind.domset}{The vertices (i.e., data points) in the dominating set of size \code{k} if it exists,
#' otherwise it yields \code{NULL}.}
#'
#' @seealso \code{\link{IndNCSdomUBtri}}, \code{\link{IndNCSdomUBTe}}, \code{\link{IndUBdom}},
#' and \code{\link{dom.exact}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#'
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$gen.points
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' IndASdomUBtri(Xp,1,Tr)
#'
#' for (k in 1:n)
#'   print(c(k,IndASdomUBtri(Xp,k,Tr,M)))
#'
#' IndASdomUBtri(Xp,k=4,Tr,M)
#'
#' P<-c(.4,.2)
#' IndASdomUBtri(P,1,Tr,M)
#'
#' IndASdomUBtri(rbind(Xp,Xp),k=2,Tr,M)
#' }
#'
#' @export IndASdomUBtri
IndASdomUBtri <- function(Xp,k,tri,M="CC")
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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  n<-nrow(Xp);
  xc<-combinat::combn(1:n,k); N1<-choose(n,k);
  xc<-matrix(xc,ncol=N1)
  dom<-0; j<-1; domset<-c();
  while (j<=N1 && dom==0)
  {
    if (IndNAStri.domset(Xp[xc[,j],],Xp,tri,M)==1)  #this is where triangle tri is used
    {dom<-1; domset<-Xp[xc[,j],];}
    j<-j+1;
  }

  list(domUB=k, #upper bound for the domination number of AS-PCD
       IndUBdom=dom, #indicator that domination number <=k
       domset=domset #a dominating set of size k (if exists)
  )
} #end of the function
#'

#################################################################

#' @title Number of arcs of Arc Slice Proximity Catch Digraphs (AS-PCDs) - multiple triangle case
#'
#' @description Returns the number of arcs and various other quantities, vectors, and lists for Arc Slice Proximity Catch Digraph
#' (AS-PCD) whose vertices are the data points in \code{Xp} in the multiple triangle case.
#'
#' AS proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points and vertex regions in each triangle are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.
#'
#' Convex hull of \code{Yp} is partitioned by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles whose union constitutes
#' the convex hull of \code{Yp} points).
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation
#' and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle or 3D point in barycentric
#' coordinates which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"} i.e., the circumcenter of each triangle.
#'
#' @return A \code{list} with the elements
#' \item{num.arcs}{Total number of arcs in all triangles, i.e., the number of arcs for the entire AS-PCD}
#' \item{num.in.conhull}{Number of \code{Xp} points in the convex hull of \code{Yp} points}
#' \item{num.in.tris}{The vector of number of \code{Xp} points in the Delaunay triangles based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the areas of Delaunay triangles based on \code{Yp} points}
#' \item{tri.num.arcs}{The \code{vector} of the number of arcs of the component of the AS-PCD in the
#' Delaunay triangles based on \code{Yp} points}
#' \item{del.tri.ind}{A matrix of indices of Delaunay triangles based on \code{Yp} points,
#' each column corresponds to the vector of indices of the vertices of one of the Delaunay triangle.}
#' \item{data.tri.ind}{A \code{vector} of indices of vertices of the Delaunay triangles in which data points reside,
#' i.e., column number of \code{del.tri.ind} for each \code{Xp} point.}
#'
#' @seealso \code{\link{NumArcsAStri}}, \code{\link{NumArcsPE}}, and \code{\link{NumArcsCS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-15; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx),runif(nx))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' plotDeltri(Xp,Yp,xlab="",ylab="")
#'
#' M<-"CC"  #try also M<-c(1,1,1)
#'
#' NumArcsAS(Xp,Yp,M)
#' NumArcsAS(Xp,Yp[1:3,],M)
#' NumArcsAS(c(.4,.2),Yp,M)
#'
#' Xp2=Xp+10
#' NumArcsAS(Xp2,Yp,M)
#'
#' NumArcsAS(Xp,rbind(Yp,Yp),M)
#' }
#'
#' @export NumArcsAS
NumArcsAS <- function(Xp,Yp,M="CC")
{
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

  if ((!is.point(M,3) && M!="CC") || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates or "CC" for circumcenter')}

  nx<-nrow(Xp); ny<-nrow(Yp)

  #Delaunay triangulation of Yp points
  Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
  Ytri<-matrix(interp::triangles(Ytrimesh)[,1:3],ncol=3); #the indices of the vertices of the Delaunay triangles (row-wise)
  ndt<-nrow(Ytri)  #number of Delaunay triangles

  inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
  NinCH<-sum(inCH)  #number of points in the convex hull

  Wvec=vector()
  for (i in 1:ndt)
  {
    ifelse(ndt==1,Tri<-Yp[Ytri,],Tri<-Yp[Ytri[i,],])  #vertices of ith triangle
    Wvec<-c(Wvec,area.polygon(Tri))
  }

  if (ny==3)
  {
    Tri.Ind<-indices.del.tri(Xp,Yp)  #returns 1's if the points Xp[i,]'s are inside triangle based on Yp, NA otherwise

    inTri<-which(Tri.Ind==1)
    NinTri<-length(inTri)  #number of points in the triangle

    if (NinTri==0)
    {Tot.Arcs<-0;
    ni.vec<-arcs<-rep(0,ndt)
    data.tri.ind =NULL
    } else
    {
      Xdt<-matrix(Xp[inTri,],ncol=2)
      tri<-as.bas.tri(Yp)$tri #convert the triangle Yp into an nonscaled basic triangle, see as.bas.tri help page
      Wvec<-area.polygon(tri)
      Tot.Arcs<-NumArcsAStri(Xdt,tri,M)$num.arcs  #number of arcs in the triangle Yp
    }
    res<-list(num.arcs=Tot.Arcs,
              tri.num.arcs=Tot.Arcs,
              num.in.conv.hull=NinTri,
              num.in.tris=NinTri,
              weight.vec=Wvec,
              del.tri.ind=t(Ytri),
              data.tri.ind=Tri.Ind )
  } else
  {
    if (NinCH==0)
    {Tot.Arcs<-0;
    ni.vec<-arcs<-rep(0,ndt)
    data.tri.ind =NULL
    } else
    {
      Tri.Ind<-indices.del.tri(Xp,Yp,Ytrimesh) #indices of triangles in which the points in the data fall

      #calculation of the total number of arcs
      ni.vec<-arcs<-vector()
      # data.del.tris = del.tris=list()
      data.tri.ind = rep(NA,nx)
      for (i in 1:ndt)
      {
        dt.ind=which(Tri.Ind==i) #which indices of data points residing in ith Delaunay triangle
        Xpi<-Xp[dt.ind,] #points in ith Delaunay triangle
        data.tri.ind[dt.ind] =i #assigning the index of the Delaunay triangle that contains the data point
        #  data.del.tris=c(data.del.tris,list(Xpi))
        ifelse(ndt==1,Tri<-Yp[Ytri,],Tri<-Yp[Ytri[i,],])  #vertices of ith triangle
        tri<-as.bas.tri(Tri)$tri #convert the triangle Tri into an nonscaled basic triangle, see as.bas.tri help page
        #  del.tris=rbind(del.tris,tri)
        ni.vec<-c(ni.vec,length(Xpi)/2)  #number of points in ith Delaunay triangle

        ifelse(identical(M,"CC"),cent<-circ.cent.tri(tri),cent<-M)
        num.arcs<-NumArcsAStri(Xpi,tri,M)$num.arcs  #number of arcs in ith triangle
        arcs<-c(arcs,num.arcs)  #number of arcs in all triangles as a vector

      }

      Tot.Arcs<-sum(arcs)  #the total number of arcs in all triangles
    }

    res<-list(num.arcs=Tot.Arcs, #number of arcs for the entire PCD
              tri.num.arcs=arcs, #vector of number of arcs for the Delaunay triangles
              num.in.conv.hull=NinCH, # number of Xp points in CH of Yp points
              num.in.tris=ni.vec, # vector of number of Xp points in the Delaunay triangles
              weight.vec=Wvec, #areas of Delaunay triangles
              del.tri.ind=t(Ytri), # indices of the Delaunay triangles, each column correponds to the indices of the vertices of one triangle
              data.tri.ind=data.tri.ind) #indices of vertices of the Delaunay triangles in which data points reside, i.e., column number of del.tri for each Xp point
  }
res

} #end of the function
#'

#################################################################

#' @title Incidence matrix for Arc Slice Proximity Catch Digraphs (AS-PCDs) - one triangle case
#'
#' @description Returns the incidence matrix for the AS-PCD whose vertices are the given 2D numerical data set, \code{Xp}.
#'
#' AS proximity regions are defined with respect to the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)} and
#' vertex regions based on the center \code{M="CC"} for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri}; default is \code{M="CC"} i.e., circumcenter of \code{tri}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#'
#' @return Incidence matrix for the AS-PCD whose vertices are 2D data set, \code{Xp},
#' and AS proximity regions are defined with respect to the triangle \code{tri} and
#' vertex regions based on circumcenter.
#'
#' @seealso \code{\link{IncMatAS}}, \code{\link{IncMatPEtri}}, and \code{\link{IncMatCStri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' IM<-IncMatAStri(Xp,Tr,M)
#' IM
#'
#' dom.greedy(IM)
#' dom.exact(IM)
#' }
#'
#' @export IncMatAStri
IncMatAStri <- function(Xp,tri,M="CC")
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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  n<-nrow(Xp)
  inc.mat<-matrix(0, nrow=n, ncol=n)
  if (n>=1)
  {
    for (i in 1:n)
    {p1<-Xp[i,]
    vert<-ifelse(isTRUE(all.equal(M,circ.cent.tri(tri)))==TRUE,rv.triCC(p1,tri)$rv,rv.tri.cent(p1,tri,M)$rv)  #vertex region for pt
    for (j in ((1:n)) )
    {p2<-Xp[j,]
    inc.mat[i,j]<-IndNAStri(p1,p2,tri,M,rv=vert)
    }
    }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title Incidence matrix for Arc Slice Proximity Catch Digraphs (AS-PCDs) - multiple triangle case
#'
#' @description Returns the incidence matrix for the AS-PCD whose vertices are a given 2D numerical data set, \code{Xp},
#' in the convex hull of \code{Yp} which is partitioned by the Delaunay triangles based on \code{Yp} points.
#'
#' AS proximity regions are defined with respect to the Delaunay triangles based on \code{Yp} points and vertex
#' regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation
#' and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle or 3D point in barycentric
#' coordinates which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"} i.e., the circumcenter of each triangle.
#'
#' @return Incidence matrix for the AS-PCD whose vertices are the 2D data set, \code{Xp},
#' and AS proximity regions are defined in the Delaunay triangles based on \code{Yp} points.
#'
#' @seealso \code{\link{IncMatAStri}}, \code{\link{IncMatPE}}, and \code{\link{IncMatCS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-"CC"  #try also M<-c(1,1,1)
#'
#' IM<-IncMatAS(Xp,Yp,M)
#' IM
#' dom.greedy(IM)  #try also dom.exact(IM)  #this might take a long time for large  nx
#'
#' IM<-IncMatAS(Xp,Yp[1:3,],M)
#'
#' IncMatAS(Xp,rbind(Yp,Yp))
#' }
#'
#' @export IncMatAS
IncMatAS <- function(Xp,Yp,M="CC")
{
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

  if (nrow(Yp)==3)
  {
    inc.mat<-IncMatAStri(Xp,Yp,M)
  } else
  {
    if ((!is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates or "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    inc.mat<-matrix(0, nrow=nx, ncol=nx)

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)  #number of Delaunay triangles

    if (nx>1)
    {
      i.tr<-rep(0,nx)  #the vector of indices for the triangles that contain the Xp points
      for (i in 1:nx)
        for (j in 1:nt)
        {
          tri<-Yp[DTr[j,],]
          if (in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri )
            i.tr[i]<-j
        }

      for (i in 1:nx)
      {p1<-Xp[i,]
      if (i.tr[i]!=0)
      {
        Yi.Tri<-Yp[DTr[i.tr[i],],]
        Yi.tri<-as.bas.tri(Yi.Tri)$tri #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.bas.tri help page
        vert<-ifelse(identical(M,"CC"),rv.triCC(p1,Yi.tri)$rv,rv.tri.cent(p1,Yi.tri,M)$rv)  #vertex region for pt
        for (j in 1:nx )
        {p2<-Xp[j,]
        inc.mat[i,j]<-IndNAStri(p1,p2,Yi.tri,M,rv=vert)
        }
      }
      }
    }

    diag(inc.mat)<-1
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for Arc Slice Proximity Catch Digraphs
#' (AS-PCDs) - standard basic triangle case
#'
#' @description Returns I(\code{p} is a dominating point of the AS-PCD) where the vertices of the AS-PCD are the 2D data set \code{Xp}, that is, returns 1 if \code{p} is a dominating
#' point of AS-PCD, returns 0 otherwise. AS proximity regions are defined with respect to the standard basic triangle,
#' \eqn{T_b}, \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' Vertex regions are based on the center \code{M="CC"} for circumcenter
#' of \eqn{T_b}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of \eqn{T_b}; default is \code{M="CC"}.
#' Point, \code{p}, is in the vertex region of vertex \code{rv} (default is \code{NULL}); vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point or not of the AS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the AS-PCD.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \eqn{T_b} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \eqn{T_b}.
#' @param rv Index of the vertex whose region contains point \code{p}, \code{rv} takes the vertex labels as \eqn{1,2,3} as
#' in the row order of the vertices in \eqn{T_b}.
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point in \code{Xp} or not
#' (default is \code{FALSE}).
#'
#' @return I(\code{p} is a dominating point of the AS-PCD) where the vertices of the AS-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Gam1AStri}} and \code{\link{Gam1PEbas.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' Gam1ASbas.tri(Xp[1,],Xp,c1,c2,M)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Gam1ASbas.tri(Xp[i,],Xp,c1,c2,M))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' #or try
#' Rv<-rv.bas.triCC(Xp[1,],c1,c2)$rv
#' Gam1ASbas.tri(Xp[1,],Xp,c1,c2,M,Rv)
#'
#' Gam1ASbas.tri(c(.2,.4),Xp,c1,c2,M)
#' Gam1ASbas.tri(c(.2,.4),c(.2,.4),c1,c2,M)
#'
#' Xp2<-rbind(Xp,c(.2,.4))
#' Gam1ASbas.tri(Xp[1,],Xp2,c1,c2,M)
#'
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tb)}
#' #need to run this when M is given in barycentric coordinates
#'
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.bas.tri(c1,c2,M)
#' }
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp)
#' points(rbind(Xp[ind.gam1,]),pch=4,col=2)
#'
#' txt<-rbind(Tb,cent,Ds)
#' xc<-txt[,1]+c(-.03,.03,.02,.06,.06,-0.05,.01)
#' yc<-txt[,2]+c(.02,.02,.03,.0,.03,.03,-.03)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' Gam1ASbas.tri(c(.4,.2),Xp,c1,c2,M)
#'
#' Gam1ASbas.tri(c(.5,.11),Xp,c1,c2,M)
#'
#' Gam1ASbas.tri(c(.5,.11),Xp,c1,c2,M,ch.data.pnt=FALSE)
#' #gives an error message if ch.data.pnt=TRUE since the point is not in the standard basic triangle
#' }
#'
#' @export
Gam1ASbas.tri <- function(p,Xp,c1,c2,M="CC",rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p))
  {stop('p must be a numeric point of dimension 2')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (isTRUE(all.equal(matrix(p,ncol=2),Xp)))
  {dom<-1; return(dom); stop}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('p is not a data point in Xp')}
  }

  A<-c(0,0); B<-c(1,0); C<-c(c1,c2); Tb<-rbind(A,B,C)
  if (in.triangle(p,Tb,boundary = TRUE)$in.tri==FALSE)
  {dom<-0; return(dom); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.bas.triCC(p,c1,c2)$rv,rv.bas.tri.cent(p,c1,c2,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  n<-nrow(Xp)
  dom<-1; i<-1;

  while (i <= n & dom==1)
  {if (IndNASbas.tri(p,Xp[i,],c1,c2,M,rv)==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for Arc Slice Proximity Catch Digraphs
#' (AS-PCDs) - one triangle case
#'
#' @description Returns I(\code{p} is a dominating point of the AS-PCD whose vertices are the 2D data set \code{Xp}),
#' that is, returns 1 if \code{p} is a dominating point of AS-PCD, returns 0 otherwise.
#' Point, \code{p}, is in the region of vertex \code{rv} (default is \code{NULL}); vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise in \code{tri}.
#'
#' AS proximity regions are defined with respect to the
#' triangle \code{tri} and vertex regions are based on the center \code{M="CC"} for circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \code{M="CC"} the circumcenter of \code{tri}.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point or not of the AS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in  barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param rv Index of the vertex whose region contains point \code{p}, \code{rv} takes the vertex labels as \eqn{1,2,3} as
#' in the row order of the vertices in \code{tri}.
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point in \code{Xp} or not
#' (default is \code{FALSE}).
#'
#' @return I(\code{p} is a dominating point of the AS-PCD whose vertices are the 2D data set \code{Xp}),
#' that is, returns 1 if \code{p} is a dominating point of the AS-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Gam1ASbas.tri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' Gam1AStri(Xp[1,],Xp,Tr,M)
#' Gam1AStri(Xp[1,],Xp[1,],Tr,M)
#' Gam1AStri(c(1.5,1.5),c(1.6,1),Tr,M)
#' Gam1AStri(c(1.6,1),c(1.5,1.5),Tr,M)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Gam1AStri(Xp[i,],Xp,Tr,M))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' #or try
#' Rv<-rv.triCC(Xp[1,],Tr)$rv
#' Gam1AStri(Xp[1,],Xp,Tr,M,Rv)
#'
#' Gam1AStri(c(.2,.4),Xp,Tr,M)
#' Gam1AStri(c(.2,.4),c(.2,.4),Tr,M)
#'
#' Xp2<-rbind(Xp,c(.2,.4))
#' Gam1AStri(Xp[1,],Xp2,Tr,M)
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#'
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.tri(Tr,M)
#' }
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp)
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(rbind(Xp[ind.gam1,]),pch=4,col=2)
#'
#' txt<-rbind(Tr,cent,Ds)
#' xc<-txt[,1]
#' yc<-txt[,2]
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' Gam1AStri(c(1.5,1.1),Xp,Tr,M)
#'
#' Gam1AStri(c(1.5,1.1),Xp,Tr,M)
#'
#' Gam1AStri(c(1.5,1.1),Xp,Tr,M,ch.data.pnt=FALSE)
#' #gives an error message if ch.data.pnt=TRUE since point p is not a data point in Xp
#' }
#'
#' @export
Gam1AStri <- function(p,Xp,tri,M="CC",rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p))
  {stop('p must be a numeric point of dimension 2')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('p must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (isTRUE(all.equal(matrix(p,ncol=2),Xp)))
  {dom<-1; return(dom); stop}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('p is not a data point in Xp')}
  }

  if (in.triangle(p,tri,boundary = TRUE)$in.tri==FALSE)
  {dom<-0; return(dom); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p,tri)$rv,rv.tri.cent(p,tri,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (IndNAStri(p,Xp[i,],tri,M,rv)==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title Check a point belong to a given data set
#'
#' @description returns \code{TRUE} if the point \code{p} of any dimension is inside the data set \code{Xp} of the same dimension as \code{p};
#' otherwise returns \code{FALSE}.
#'
#' @param p A 2D point for which the function checks membership to the data set \code{Xp}.
#' @param Xp A set of 2D points representing the set of data points.
#'
#' @return \code{TRUE} if \code{p} belongs to the data set \code{Xp}.
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-10
#' Xp<-cbind(runif(n),runif(n))
#'
#' P<-Xp[7,]
#' is.in.data(P,Xp)
#' is.in.data(P,Xp[7,])
#'
#' P<-Xp[7,]+10^(-7)
#' is.in.data(P,Xp)
#'
#' P<-Xp[7,]+10^(-9)
#' is.in.data(P,Xp)
#'
#' is.in.data(P,P)
#'
#' is.in.data(c(2,2),c(2,2))
#'
#' #for 1D data
#' n<-10
#' Xp<-runif(n)
#'
#' P<-Xp[7]
#' is.in.data(P,Xp[7])  #this works because both entries are treated as 1D vectors but
#' #is.in.data(P,Xp) does not work since entries are treated as vectors of different dimensions
#'
#' Xp<-as.matrix(Xp)
#' is.in.data(P,Xp)
#' #this works, because P is a 1D point, and Xp is treated as a set of 10 1D points
#'
#' P<-Xp[7]+10^(-7)
#' is.in.data(P,Xp)
#'
#' P<-Xp[7]+10^(-9)
#' is.in.data(P,Xp)
#'
#' is.in.data(P,P)
#'
#' #for 3D data
#' n<-10
#' Xp<-cbind(runif(n),runif(n),runif(n))
#'
#' P<-Xp[7,]
#' is.in.data(P,Xp)
#' is.in.data(P,Xp[7,])
#'
#' P<-Xp[7,]+10^(-7)
#' is.in.data(P,Xp)
#'
#' P<-Xp[7,]+10^(-9)
#' is.in.data(P,Xp)
#'
#' is.in.data(P,P)
#'
#' n<-10
#' Xp<-cbind(runif(n),runif(n))
#' P<-Xp[7,]
#' is.in.data(P,Xp)
#' }
#'
#' @export is.in.data
is.in.data <- function(p,Xp)
{
  if (!is.numeric(as.matrix(p)) || !is.numeric(as.matrix(Xp)) )
  {stop ('p and Xp must be numeric')}

  if (!is.vector(p))
  {stop('p must be a vector')}

  dimp<-length(p)

  ins<-FALSE
  if (is.vector(Xp))
  {dimXp<-length(Xp);
  if (dimp != dimXp )
  {stop('Both arguments must be of the same dimension')
  } else
  {
    if (isTRUE(all.equal(p,Xp)))
      ins<-TRUE }
  } else
  {
    if (!is.matrix(Xp) && !is.data.frame(Xp))
    {stop('Xp must be a matrix or a data frame, each row representing a point')}
    Xp<-as.matrix(Xp)
    dimXp<-ncol(Xp);
    if (dimp != dimXp)
    {stop('p and Xp must be of the same dimension')}
    nXp<-nrow(Xp)
    i<-1
    cnt<-0

    while (i <= nXp & cnt==0)
    {
      if (isTRUE(all.equal(p,Xp[i,])))
      {ins<-TRUE; cnt<-1}
      else
      {i<-i+1}
    }
  }
  ins
} #end of the function
#'

#################################################################

#' @title The indicator for two points being a dominating set for Arc Slice Proximity Catch Digraphs
#' (AS-PCDs) - standard basic triangle case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is a dominating set of AS-PCD\eqn{)} where vertices of AS-PCD are the 2D
#' data set \code{Xp}), that is, returns 1 if \{\code{p1,p2}\} is a dominating set of AS-PCD, returns 0 otherwise.
#'
#' AS proximity regions are defined with respect to the standard basic triangle \eqn{T_b=T(c(0,0),c(1,0),c(c1,c2))},
#' In the standard basic triangle, \eqn{T_b}, \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' Point, \code{p1}, is in the vertex region of vertex \code{rv1} (default is \code{NULL})
#' and point, \code{p2}, is in the vertex region of vertex \code{rv2} (default is \code{NULL}); vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise.
#'
#' Vertex regions are based on the center \code{M="CC"} for circumcenter
#' of \eqn{T_b}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of \eqn{T_b}; default is \code{M="CC"}.
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1} and \code{p2} are data points in \code{Xp} or not
#' (default is \code{FALSE}), so by default this function checks whether the points \code{p1} and \code{p2} would be a
#' dominating set if they actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p1,p2 Two 2D points to be tested for constituting a dominating set of the AS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the AS-PCD.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \eqn{T_b} or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \eqn{T_b}.
#' @param rv1,rv2 The indices of the vertices whose regions contains \code{p1} and \code{p2}, respectively.
#' They take the vertex labels as \eqn{1,2,3} as in the row order of the vertices in \eqn{T_b}
#' (default is \code{NULL} for both).
#' @param ch.data.pnts A logical argument for checking whether points \code{p1} and \code{p2} are data points in \code{Xp} or not
#' (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the AS-PCD\eqn{)} where the vertices of AS-PCD are the 2D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of AS-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Gam2AStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' Gam2ASbas.tri(Xp[1,],Xp[2,],Xp,c1,c2,M)
#' Gam2ASbas.tri(Xp[1,],Xp[1,],Xp,c1,c2,M)  #one point can not a dominating set of size two
#'
#' Gam2ASbas.tri(c(.2,.4),c(.2,.5),rbind(c(.2,.4),c(.2,.5)),c1,c2,M)
#'
#' ind.gam2<-vector()
#' for (i in 1:(n-1))
#'   for (j in (i+1):n)
#'   {if (Gam2ASbas.tri(Xp[i,],Xp[j,],Xp,c1,c2,M)==1)
#'    ind.gam2<-rbind(ind.gam2,c(i,j))}
#' ind.gam2
#'
#' #or try
#' rv1<-rv.bas.triCC(Xp[1,],c1,c2)$rv
#' rv2<-rv.bas.triCC(Xp[2,],c1,c2)$rv
#' Gam2ASbas.tri(Xp[1,],Xp[2,],Xp,c1,c2,M,rv1,rv2)
#' Gam2ASbas.tri(c(.2,.4),Xp[2,],Xp,c1,c2,M,rv1,rv2)
#'
#' #or try
#' rv1<-rv.bas.triCC(Xp[1,],c1,c2)$rv
#' Gam2ASbas.tri(Xp[1,],Xp[2,],Xp,c1,c2,M,rv1)
#'
#' #or try
#' Rv2<-rv.bas.triCC(Xp[2,],c1,c2)$rv
#' Gam2ASbas.tri(Xp[1,],Xp[2,],Xp,c1,c2,M,rv2=Rv2)
#'
#' Gam2ASbas.tri(c(.3,.2),c(.35,.25),Xp,c1,c2,M)
#' }
#'
#' @export
Gam2ASbas.tri <- function(p1,p2,Xp,c1,c2,M="CC",rv1=NULL,rv2=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (isTRUE(all.equal(matrix(rbind(p1,p2),ncol=2),Xp)))
  {dom<-1; return(dom); stop}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  A<-c(0,0); B<-c(1,0); C<-c(c1,c2); Tb<-rbind(A,B,C);

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points, p1 and p2, are data points in Xp')}
  }

  if (isTRUE(all.equal(p1,p2)))
  {dom<-0; return(dom); stop}

  if (is.null(rv1))
  { rv1<-ifelse(isTRUE(all.equal(M,CC)),rv.bas.triCC(p1,c1,c2)$rv,rv.bas.tri.cent(p1,c1,c2,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv1) || sum(rv1==c(1,2,3))!=1)
  {stop('vertex index, rv1, must be 1, 2 or 3')}}

  if (is.null(rv2))
  { rv2<-ifelse(isTRUE(all.equal(M,circ.cent.tri(Tb)))== TRUE,rv.bas.triCC(p2,c1,c2)$rv,rv.bas.tri.cent(p2,c1,c2,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv2) || sum(rv2==c(1,2,3))!=1)
  {stop('vertex index, rv2, must be 1, 2 or 3')}}

  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IndNASbas.tri(p1,Xp[i,],c1,c2,M,rv=rv1),IndNASbas.tri(p2,Xp[i,],c1,c2,M,rv=rv2))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for two points constituting a dominating set for Arc Slice Proximity Catch Digraphs
#' (AS-PCDs) - one triangle case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is a dominating set of the AS-PCD\eqn{)} where vertices of the AS-PCD are
#' the 2D data set \code{Xp}), that is, returns 1 if \{\code{p1,p2}\} is a dominating set of AS-PCD, returns 0 otherwise.
#'
#' AS proximity regions are defined with respect to the triangle \code{tri}.
#' Point, \code{p1}, is in the region of vertex \code{rv1} (default is \code{NULL})
#' and point, \code{p2}, is in the region of vertex \code{rv2} (default is \code{NULL}); vertices (and hence \code{rv1} and \code{rv2})
#' are labeled as \eqn{1,2,3} in the order they are stacked row-wise in \code{tri}.
#'
#' Vertex regions are based on
#' the center \code{M="CC"} for circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \code{M="CC"} the circumcenter of \code{tri}.
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1} and \code{p2} are data
#' points in \code{Xp} or not (default is \code{FALSE}), so by default this function checks whether the points \code{p1} and \code{p2}
#' would constitute dominating set if they actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p1,p2 Two 2D points to be tested for constituting a dominating set of the AS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or
#' a 3D point in  barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param rv1,rv2 The indices of the vertices whose regions contains \code{p1} and \code{p2}, respectively.
#' They take the vertex labels as \eqn{1,2,3} as in the row order of the vertices in \code{tri}
#' (default is \code{NULL} for both).
#' @param ch.data.pnts A logical argument for checking whether points \code{p1} and \code{p2} are data points in \code{Xp} or not
#' (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the AS-PCD\eqn{)} where vertices of the AS-PCD are the 2D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of AS-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Gam2ASbas.tri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' Gam2AStri(Xp[1,],Xp[2,],Xp,Tr,M)
#' Gam2AStri(Xp[1,],Xp[1,],Xp,Tr,M)  #same two points cannot be a dominating set of size 2
#'
#' Gam2AStri(c(.2,.4),Xp[2,],Xp,Tr,M)
#' Gam2AStri(c(.2,.4),c(.2,.5),Xp,Tr,M)
#' Gam2AStri(c(.2,.4),c(.2,.5),rbind(c(.2,.4),c(.2,.5)),Tr,M)
#'
#' #or try
#' rv1<-rv.triCC(c(.2,.4),Tr)$rv
#' rv2<-rv.triCC(c(.2,.5),Tr)$rv
#' Gam2AStri(c(.2,.4),c(.2,.5),rbind(c(.2,.4),c(.2,.5)),Tr,M,rv1,rv2)
#'
#' ind.gam2<-vector()
#' for (i in 1:(n-1))
#'   for (j in (i+1):n)
#'   {if (Gam2AStri(Xp[i,],Xp[j,],Xp,Tr,M)==1)
#'    ind.gam2<-rbind(ind.gam2,c(i,j))}
#' ind.gam2
#'
#' #or try
#' rv1<-rv.triCC(Xp[1,],Tr)$rv
#' rv2<-rv.triCC(Xp[2,],Tr)$rv
#' Gam2AStri(Xp[1,],Xp[2,],Xp,Tr,M,rv1,rv2)
#'
#' #or try
#' rv1<-rv.triCC(Xp[1,],Tr)$rv
#' Gam2AStri(Xp[1,],Xp[2,],Xp,Tr,M,rv1)
#'
#' #or try
#' Rv2<-rv.triCC(Xp[2,],Tr)$rv
#' Gam2AStri(Xp[1,],Xp[2,],Xp,Tr,M,rv2=Rv2)
#'
#' Gam2AStri(c(1.3,1.2),c(1.35,1.25),Xp,Tr,M)
#' }
#'
#' @export Gam2AStri
Gam2AStri <- function(p1,p2,Xp,tri,M="CC",rv1=NULL,rv2=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (isTRUE(all.equal(matrix(rbind(p1,p2),ncol=2),Xp)))
  {dom<-1; return(dom); stop}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points, p1 and p2, are data points in Xp')}
  }

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {dom<-0; return(dom); stop}

  if (is.null(rv1))
  { rv1<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p1,tri)$rv,rv.tri.cent(p1,tri,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv1) || sum(rv1==c(1,2,3))!=1)
  {stop('vertex index, rv1, must be 1, 2 or 3')}}

  if (is.null(rv2))
  { rv2<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p2,tri)$rv,rv.tri.cent(p2,tri,M)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(rv2) || sum(rv2==c(1,2,3))!=1)
  {stop('vertex index, rv2, must be 1, 2 or 3')}}

  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IndNAStri(p1,Xp[i,],tri,M,rv=rv1),IndNAStri(p2,Xp[i,],tri,M,rv=rv2))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The arcs of Arc Slice Proximity Catch Digraph (AS-PCD) for 2D data - one triangle case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for data set \code{Xp} as the vertices
#' of AS-PCD and related parameters and the quantities of the digraph.
#'
#' AS proximity regions are constructed with respect to the triangle \code{tri}, i.e.,
#' arcs may exist for points only inside \code{tri}.
#' It also provides various descriptions and quantities about the arcs of the AS-PCD
#' such as number of arcs, arc density, etc.
#'
#' Vertex regions are based on the center \code{M="CC"} for
#' circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric
#' coordinates in the interior of the triangle \code{tri}; default is \code{M="CC"} the circumcenter of \code{tri}.
#' The different consideration of circumcenter vs any other interior center of the triangle is because
#' the projections from circumcenter are orthogonal to the edges, while projections of \code{M} on the edges are on the extensions
#' of the lines connecting \code{M} and the vertices.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, it is the center used to construct the vertex regions.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support triangle.}
#' \item{tess.name}{Name of data set used in tessellation (i.e., vertices of the triangle).}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of AS-PCD for 2D data set \code{Xp} as vertices of the digraph}
#' \item{E}{Heads (or arrow ends) of the arcs of AS-PCD for 2D data set \code{Xp} as vertices of the digraph}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{ArcsAS}}, \code{\link{ArcsPEtri}}, \code{\link{ArcsCStri}}, \code{\link{ArcsPE}},
#' and \code{\link{ArcsCS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#'
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also  M<-c(1.6,1.2) or M<-circ.cent.tri(Tr)
#'
#' Arcs<-ArcsAStri(Xp,Tr,M) #try also Arcs<-ArcsAStri(Xp,Tr)
#' #uses the default center, namely circumcenter for M
#' Arcs
#' summary(Arcs)
#' plot(Arcs) #use plot(Arcs,asp=1) if M=CC
#'
#' #can add vertex regions
#' #but we first need to determine center is the circumcenter or not,
#' #see the description for more detail.
#' CC<-circ.cent.tri(Tr)
#' M = as.numeric(Arcs$parameters[[1]])
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.tri(Tr,M)
#' }
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' #now we add the vertex names and annotation
#' txt<-rbind(Tr,cent,Ds)
#' xc<-txt[,1]+c(-.02,.03,.02,.03,.04,-.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.03,.06,.04,.05,-.07)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export ArcsAStri
ArcsAStri <- function(Xp,tri,M="CC")
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(tri))

  if (!is.numeric(as.matrix(Xp)) )
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

  CC <- circ.cent.tri(tri)
  if (identical(M,"CC"))
  {M<-CC
  } else
  { if (!is.point(M) && !is.point(M,3))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

    if (dimension(M)==3)
    {M<-bary2cart(M,tri)}

    if ( !(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri) )
    {stop('center is not the circumcenter or not in the interior of the triangle')}
  }

  n<-nrow(Xp)
  in.tri<-rep(0,n)
  for (i in 1:n)
    in.tri[i]<-in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri #indices of the Xp points inside the triangle

  Xtri<-Xp[in.tri==1,] #the Xp points inside the triangle
  n2<-length(Xtri)/2

  #the arcs of AS-PCDs
  S<-E<-NULL #S is for source and E is for end points for the arcs
  if (n2>1)
  {
    for (j in 1:n2)
    {
      p1<-Xtri[j,];
      RV1<-ifelse(isTRUE(all.equal(M,CC)),rv.triCC(p1,tri)$rv,rv.tri.cent(p1,tri,M)$rv) #vertex region for p1
      for (k in (1:n2)[-j])  #to avoid loops
      {
        p2<-Xtri[k,];
        if (IndNAStri(p1,p2,tri,M,RV1)==1)
        {
          S <-rbind(S,Xtri[j,]); E <-rbind(E,Xtri[k,]);
        }
      }
    }
  }

  param<-list(M)
  Mr<-round(M,2)
  if (identical(M,"CC") || isTRUE(all.equal(M,CC)))
      { cname <-"CC"
      names(param)<-c("circumcenter")
        typ<-"Arc Slice Proximity Catch Digraph (AS-PCD) for 2D Points in the Triangle with CC-Vertex Regions"
        main.txt<-"Arcs of AS-PCD with CC-Vertex Regions"
        } else
        {
          cname <-"M"
          names(param)<-c("center")
        typ<-paste("Arc Slice Proximity Catch Digraph (AS-PCD) for 2D Points in the Triangle with Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
        main.txt<-paste("Arcs of AS-PCD with Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
        }

  nvert<-n2; ny<-3; ntri<-1; narcs=ifelse(!is.null(S),nrow(S),0);
  arc.dens<-ifelse(nvert>1,narcs/(nvert*(nvert-1)),NA)

  quantities<-c(nvert,ny,ntri,narcs,arc.dens)
  names(quantities)<-c("number of vertices", "number of partition points",
                       "number of triangles","number of arcs", "arc density")

  res<-list(
    type=typ,
    parameters=param,
    tess.points=tri, tess.name=yname, #tessellation points
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

#' @title The plot of the arcs of Arc Slice Proximity Catch Digraph (AS-PCD) for a 2D data set - one triangle case
#'
#' @description Plots the arcs of AS-PCD whose vertices are the data points, \code{Xp} and the triangle \code{tri}. AS proximity regions
#' are constructed with respect to the triangle \code{tri}, i.e., only for \code{Xp} points inside the triangle \code{tri}.
#'
#' Vertex regions are based on the center \code{M="CC"} for
#' circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric
#' coordinates in the interior of the triangle \code{tri}; default is \code{M="CC"} the circumcenter of \code{tri}.
#' When the center is the circumcenter, \code{CC}, the vertex regions are constructed based on the
#' orthogonal projections to the edges, while with any interior center \code{M}, the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param asp A \code{numeric} value, giving the aspect ratio for \eqn{y} axis to \eqn{x}-axis \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param vert.reg A logical argument to add vertex regions to the plot, default is \code{vert.reg=FALSE}.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of the AS-PCD for a 2D data set \code{Xp} where AS proximity regions
#' are defined with respect to the triangle \code{tri}; also plots the triangle \code{tri}
#'
#' @seealso \code{\link{plotASarcs}}, \code{\link{plotPEarcs.tri}}, \code{\link{plotPEarcs}},
#' \code{\link{plotCSarcs.tri}}, and \code{\link{plotCSarcs}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g  #try also Xp<-cbind(runif(n,1,2),runif(n,0,2))
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also  #M<-c(1.6,1.2)
#'
#' plotASarcs.tri(Xp,Tr,M,main="Arcs of AS-PCD",xlab="",ylab="")
#'
#' plotASarcs.tri(Xp,Tr,M,main="Arcs of AS-PCD",xlab="",ylab="",vert.reg = TRUE)
#'
#' # or try the default center
#' #plotASarcs.tri(Xp,Tr,asp=1,main="arcs of AS-PCD",xlab="",ylab="",vert.reg = TRUE);
#' #M = (ArcsAStri(Xp,Tr)$param)$c #the part "M = as.numeric(ArcsAStri(Xp,Tr)$param)" is optional,
#' #for the below annotation of the plot
#'
#' #can add vertex labels and text to the figure (with vertex regions)
#' #but first we need to determine whether the center used for vertex regions is CC or not
#' #see the description for more detail.
#' CC<-circ.cent.tri(Tr)
#'
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.tri(Tr,M)
#' }
#'
#' #now we add the vertex names and annotation
#' txt<-rbind(Tr,cent,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.01,.05,-0.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.02,.07,.02,.05,-.06)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export plotASarcs.tri
plotASarcs.tri <- function(Xp,tri,M="CC",asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,vert.reg=FALSE,...)
{
  ArcsAS<-ArcsAStri(Xp,tri,M)
  S<-ArcsAS$S
  E<-ArcsAS$E
  cent = (ArcsAS$param)$c

  if (is.null(xlim))
  {Xlim<-range(tri[,1],Xp[,1],cent[1])
  xd<-Xlim[2]-Xlim[1]
  xlim=Xlim+xd*c(-.05,.05)
  }
  if (is.null(ylim))
  {Ylim<-range(tri[,2],Xp[,2],cent[2])
  yd<-Ylim[2]-Ylim[1]
  ylim=Ylim+yd*c(-.05,.05)
  }

  if ( isTRUE(all.equal(cent,circ.cent.tri(tri))) )
  {M="CC"}

  if (is.null(main))
  {if (identical(M,"CC")){
    main="Arcs of AS-PCD with CC-Vertex Regions"
  } else {Mr=round(cent,2);
  Mvec= paste(Mr, collapse=",")
  main=paste("Arcs of AS-PCD with M = (",Mvec,")",sep="")}
  }

  if (vert.reg)
  {main=c(main,"\n (vertex regions added)")}

  plot(Xp,main=main, asp=asp, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3,...)
  polygon(tri,...)
  if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}
  if(vert.reg){
    cent = (ArcsAS$par)$c
    CC<-circ.cent.tri(tri)
    A = tri[1,]; B = tri[2,]; C = tri[3,]

    if (isTRUE(all.equal(cent,CC)) || identical(M,"CC"))
    {D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
    Ds<-rbind(D1,D2,D3)
    } else
    { Ds<-cent2edges.tri(tri,M)}
    L<-rbind(cent,cent,cent); R<-Ds
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }
} #end of the function
#'

#################################################################

#' @title The plot of the Arc Slice (AS) Proximity Regions for a 2D data set - one triangle case
#'
#' @description Plots the points in and outside of the triangle \code{tri} and also the AS proximity regions
#' for points in data set \code{Xp}.
#'
#' AS proximity regions are defined with respect to the triangle \code{tri},
#' so AS proximity regions are defined only for points inside the triangle \code{tri} and
#' vertex regions are based on the center \code{M="CC"} for
#' circumcenter of \code{tri}; or \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric
#' coordinates in the interior of the triangle \code{tri}; default is \code{M="CC"} the circumcenter of \code{tri}.
#' When vertex regions are constructed with circumcenter, \code{CC}, the vertex regions are constructed based on the
#' orthogonal projections to the edges, while with any interior center \code{M}, the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points for which AS proximity regions are constructed.
#' @param tri Three 2D points, stacked row-wise, each row representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of the triangle \code{tri} or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"} i.e., the circumcenter of \code{tri}.
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param vert.reg A logical argument to add vertex regions to the plot, default is \code{vert.reg=FALSE}.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the AS proximity regions for points inside the triangle \code{tri} (and only the points outside \code{tri})
#'
#' @seealso \code{\link{plotASregs}}, \code{\link{plotPEregs.tri}}, \code{\link{plotPEregs}},
#' \code{\link{plotCSregs.tri}}, and \code{\link{plotCSregs}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp0<-runif.tri(n,Tr)$g
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also  #M<-c(1.6,1.2);
#'
#' plotASregs.tri(Xp0,Tr,M,main="Proximity Regions for AS-PCD", xlab="",ylab="")
#' Xp = Xp0[1,]
#' plotASregs.tri(Xp,Tr,M,main="Proximity Regions for AS-PCD", xlab="",ylab="")
#'
#' #can plot the arcs of the AS-PCD
#' #plotASarcs.tri(Xp,Tr,M,main="Arcs of AS-PCD",xlab="",ylab="")
#'
#' plotASregs.tri(Xp,Tr,M,main="Proximity Regions for AS-PCD", xlab="",ylab="",vert.reg=TRUE)
#'
#' # or try the default center
#' #plotASregs.tri(Xp,Tr,main="Proximity Regions for AS-PCD", xlab="",ylab="",vert.reg=TRUE);
#' M = (ArcsAStri(Xp,Tr)$param)$c #the part "M = as.numeric(ArcsAStri(Xp,Tr)$param)" is optional,
#' #for the below annotation of the plot
#'
#' #can add vertex labels and text to the figure (with vertex regions)
#' #but first we need to determine whether the center used for vertex regions is CC or not
#' #see the description for more detail.
#' CC<-circ.cent.tri(Tr)
#' #Arcs<-ArcsAStri(Xp,Tr,M)
#' #M = as.numeric(Arcs$parameters)
#' if (isTRUE(all.equal(M,CC)) || identical(M,"CC"))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-cent2edges.tri(Tr,M)
#' }
#'
#' #now we add the vertex names and annotation
#' txt<-rbind(Tr,cent,Ds)
#' xc<-txt[,1]+c(-.02,.03,.03,.03,.05,-0.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.02,.07,.02,.05,-.06)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export plotASregs.tri
plotASregs.tri <- function(Xp,tri,M="CC",main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,vert.reg=FALSE,...)
{
  if (!is.numeric(as.matrix(Xp)) )
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

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") ||isTRUE(all.equal(M,CC)) )
  {cent=CC; M="CC"
  } else if (dimension(M)==3)
  { M<-cent<-bary2cart(M,tri)
  } else
  {cent=M}

    if (!identical(M,"CC") && in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
    {stop('center is not the circumcenter or not in the interior of the triangle')}

  nx<-nrow(Xp)
  in.tri<-rep(0,nx)
  for (i in 1:nx)
    in.tri[i]<-in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri #indices of the Xp points inside the triangle

  Xtri<-matrix(Xp[in.tri==1,],ncol=2)  #the Xp points inside the triangle
  nx2<-nrow(Xtri)

  if (is.null(xlim))
  {Xlim<-range(tri[,1],Xp[,1],cent[1])
  xd<-Xlim[2]-Xlim[1]
  xlim=Xlim+xd*c(-.05,.05)
  }
  if (is.null(ylim))
  {Ylim<-range(tri[,2],Xp[,2],cent[2])
  yd<-Ylim[2]-Ylim[1]
  ylim=Ylim+yd*c(-.05,.05)
  }

  if (is.null(main))
  {if (identical(M,"CC")) {
    main="AS Proximity Regions with CC-Vertex Regions"
  } else {  Mr=round(M,2)
  Mvec= paste(Mr, collapse=",")
  main=paste("AS Proximity Regions with M = (",Mvec,")",sep="")}
  }

  if (vert.reg)
  {main=c(main,"\n (vertex regions added)")}

  plot(Xp,main=main, asp=1, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3,...)
  polygon(tri,lty=2)
  if (nx2>=1)
  {
    for (i in 1:nx2)
    {
      P1<-Xtri[i,]
      rv<-ifelse(identical(M,"CC"),rv.triCC(P1,tri)$rv,rv.tri.cent(P1,tri,M)$rv)  #vertex region for pt
      RV<-tri[rv,]
      rad<-Dist(P1,RV)

      Int.Pts<-NAStri(P1,tri,M)
      L<-matrix(Int.Pts$L,ncol=2); R<-matrix(Int.Pts$R,ncol=2)
      segments(L[,1], L[,2], R[,1], R[,2], lty=1,col=2)
      Arcs<-Int.Pts$a;
      if (!is.null(Arcs))
      {
        K<-nrow(Arcs)/2
        for (j in 1:K)
        {A1<-Arcs[2*j-1,]; A2<-Arcs[2*j,];
        angles<-angle.str2end(A1,P1,A2)$c

        test.ang1<-angles[1]+(.01)*(angles[2]-angles[1])
        test.Pnt<-P1+rad*c(cos(test.ang1),sin(test.ang1))
        if (!in.triangle(test.Pnt,tri,boundary = TRUE)$i) {angles<-c(min(angles),max(angles)-2*pi)}
        plotrix::draw.arc(P1[1],P1[2],rad,angle1=angles[1],angle2=angles[2],col=2)
        }
      }
    }
  }

  if(vert.reg){
    A = tri[1,]; B = tri[2,]; C = tri[3,]
    if (identical(M,"CC"))
    {D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
    Ds<-rbind(D1,D2,D3)
    } else
    {
      cent=M
      Ds<-cent2edges.tri(tri,M)}
    L<-rbind(cent,cent,cent); R<-Ds
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }

} #end of the function
#'

#################################################################

#' @title The arcs of Arc Slice Proximity Catch Digraph (AS-PCD) for a 2D data set - multiple triangle case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) of AS-PCD whose vertices are the data set \code{Xp}
#' and related parameters and the quantities of the digraph.
#'
#' AS proximity regions are defined with respect to the Delaunay triangles based on
#' \code{Yp} points, i.e., AS proximity regions are defined only for \code{Xp} points
#' inside the convex hull of \code{Yp} points.
#' That is, arcs may exist for points only inside the convex hull of \code{Yp} points.
#' It also provides various descriptions and quantities about the arcs of the AS-PCD
#' such as number of arcs, arc density, etc.
#'
#' Vertex regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.
#' \code{M} must be entered in barycentric coordinates unless it is the circumcenter.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangulation. The Delaunay
#' triangles partition the convex hull of \code{Yp} points.
#' @param M The center of the triangle. \code{"CC"} represents the circumcenter of each Delaunay triangle
#' or 3D point in barycentric coordinates which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"} i.e., the circumcenter of each triangle. \code{M} must be entered in barycentric coordinates
#' unless it is the circumcenter.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, it is the center used to construct the vertex regions,
#' default is circumcenter, denoted as \code{"CC"}, otherwise given in barycentric coordinates.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the Delaunay triangulation based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation, i.e., \code{Yp}}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of AS-PCD for 2D data set \code{Xp} in the multiple triangle case
#' as the vertices of the digraph}
#' \item{E}{Heads (or arrow ends) of the arcs of AS-PCD for 2D data set \code{Xp} in the multiple triangle case
#' as the vertices of the digraph}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of intervals, number of arcs, and arc density.}
#'
#' @seealso \code{\link{ArcsAStri}}, \code{\link{ArcsPEtri}}, \code{\link{ArcsCStri}}, \code{\link{ArcsPE}},
#' and \code{\link{ArcsCS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-5;  #try also nx=20; nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' Arcs<-ArcsAS(Xp,Yp,M) #try also the default M with Arcs<-ArcsAS(Xp,Yp)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' ArcsAS(Xp,Yp[1:3,],M)
#' }
#'
#' @export ArcsAS
ArcsAS <- function(Xp,Yp,M="CC")
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
  {stop('Xp must be of dimension n x 2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (nrow(Yp)==3)
  {
    res<-ArcsAStri(Xp,Yp,M)
  } else {
    if ((!is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates or "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)  #the Xp points inside the convex hull of Yp points

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)
    nx2<-nrow(Xch)  #number of Xp points inside the convex hull of Yp points

    S<-E<-NULL #S is for source and E is for end points for the arcs
    if (nx2>1)
    {
      i.tr<-rep(0,nx2)  #the vector of indices for the triangles that contain the Xp points
      for (i in 1:nx2)
        for (j in 1:nt)
        {
          tri<-Yp[DTr[j,],]
          if (in.triangle(Xch[i,],tri,boundary=TRUE)$in.tri )
            i.tr[i]<-j
        }

      for (i in 1:nt)
      {
        Xl<-matrix(Xch[i.tr==i,],ncol=2)
        if (nrow(Xl)>1)
        {
          Yi.Tri<-Yp[DTr[i,],] #vertices of the ith triangle
          Yi.tri<-as.bas.tri(Yi.Tri)$tri #convert the triangle Tri into an nonscaled basic triangle, see as.bas.tri help page
          nl<-nrow(Xl)
          ifelse(identical(M,"CC"), rv.ind<-rel.verts.triCC(Xl,tri=Yi.tri)$rv,
                 rv.ind<-rel.verts.tri.cent(Xl,tri=Yi.tri,M)$rv)  #vertex region for pt
          for (j in 1:nl)
          {RV<-rv.ind[j]
          for (k in (1:nl)[-j])  # to avoid loops
            if (IndNAStri(Xl[j,],Xl[k,],Yi.tri,M,rv=RV)==1 )
            {
              S <-rbind(S,Xl[j,]); E <-rbind(E,Xl[k,]);
            }
          }
        }
      }
    }

   cname <-ifelse(identical(M,"CC"),"CC","M")

    param<-list(M)
    names(param)<-c("center")

    Mvec <- paste(M, collapse=",");
    if (identical(M,"CC"))
        {
          typ<-"Arc Slice Proximity Catch Digraph (AS-PCD) for 2D Points in Multiple Triangles with CC-Vertex Regions"
          main.txt<-"Arcs of AS-PCD with CC-Vertex Regions"
          } else
          { typ<-paste("Arc Slice Proximity Catch Digraph (AS-PCD) for 2D Points in Multiple Triangles with Center ", cname," = (",Mvec,")",sep="")
      main.txt <- paste("Arcs of AS-PCD with Center ", cname," = (",Mvec,")",sep="")
    }

    nvert<-nx2; ny<-nrow(Yp); ntri<-nt; narcs=ifelse(!is.null(S),nrow(S),0);
    arc.dens<-ifelse(nvert>1,narcs/(nvert*(nvert-1)),NA)

    quantities<-c(nvert,ny,ntri,narcs,arc.dens)
    names(quantities)<-c("number of vertices", "number of partition points",
                         "number of triangles","number of arcs", "arc density")

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
  }
  res
} #end of the function
#'

#################################################################

#' @title The plot of the arcs of Arc Slice Proximity Catch Digraph (AS-PCD) for a 2D data set -
#' multiple triangle case
#'
#' @description Plots the arcs of AS-PCD whose vertices are the data points in \code{Xp} and Delaunay triangles based on \code{Yp} points.
#'
#' AS proximity regions are constructed with respect to the Delaunay triangles based on \code{Yp} points, i.e.,
#' AS proximity regions are defined only for \code{Xp} points inside the convex hull of \code{Yp} points.
#' That is, arcs may exist for \code{Xp} points only inside the convex hull of \code{Yp} points.
#'
#' Vertex regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the AS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangulation. The Delaunay
#' triangles partition the convex hull of \code{Yp} points.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle
#' or 3D point in barycentric coordinates which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"} i.e., the circumcenter of each triangle.
#' @param asp A \code{numeric} value, giving the aspect ratio for \eqn{y} axis to \eqn{x}-axis \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of the AS-PCD for a 2D data set \code{Xp} where AS proximity regions
#' are defined with respect to the Delaunay triangles based on \code{Yp} points; also plots the Delaunay triangles
#' based on \code{Yp} points.
#'
#' @seealso \code{\link{plotASarcs.tri}}, \code{\link{plotPEarcs.tri}}, \code{\link{plotPEarcs}},
#' \code{\link{plotCSarcs.tri}}, and \code{\link{plotCSarcs}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' #plotASarcs(Xp,Yp,M,xlab="",ylab="")
#' plotASarcs(Xp,Yp,M,asp=1,xlab="",ylab="")
#'
#' plotASarcs(Xp,Yp[1:3,],M,xlab="",ylab="")
#' }
#'
#' @export plotASarcs
plotASarcs <- function(Xp,Yp,M="CC",asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
{
  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (nrow(Yp)==3)
  {
    plotASarcs.tri(Xp,Yp,M,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    ArcsAS<-ArcsAS(Xp,Yp,M)
    S<-ArcsAS$S
    E<-ArcsAS$E

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    if (is.null(xlim))
    {xlim<-range(Yp[,1],Xp[,1])}
    if (is.null(ylim))
    {ylim<-range(Yp[,2],Xp[,2])}

    xr<-xlim[2]-xlim[1]
    yr<-ylim[2]-ylim[1]

    if (is.null(main))
    {if (identical(M,"CC")){
      main="Arcs of AS-PCD with CC-Vertex Regions"
    } else {Mvec= paste(M, collapse=",")
    main=paste("Arcs of AS-PCD with M = (",Mvec,")",sep="")}
    }

    plot(rbind(Xp),asp=asp,main=main, xlab=xlab, ylab=ylab,xlim=xlim+xr*c(-.05,.05),
         ylim=ylim+yr*c(-.05,.05),pch=".",cex=3,...)
    interp::plot.triSht(DTmesh, add=TRUE, do.points = TRUE,...)
    if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}
  }
} #end of the function
#'

#################################################################

#' @title The plot of the Arc Slice (AS) Proximity Regions for a 2D data set - multiple triangle case
#'
#' @description Plots the \code{Xp} points in and outside of the convex hull of \code{Yp} points and also plots the AS proximity regions
#' for \code{Xp} points and Delaunay triangles based on \code{Yp} points.
#'
#' AS proximity regions are constructed with respect
#' to the Delaunay triangles based on \code{Yp} points (these triangles partition the convex hull of \code{Yp} points),
#' i.e., AS proximity regions are only defined for \code{Xp} points inside the convex hull of \code{Yp} points.
#'
#' Vertex regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle; default is \code{M="CC"} i.e., circumcenter of each triangle.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}) for more on AS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points for which AS proximity regions are constructed.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangulation. The Delaunay
#' triangles partition the convex hull of \code{Yp} points.
#' @param M The center of the triangle. \code{"CC"} stands for circumcenter of each Delaunay triangle or 3D point in barycentric
#' coordinates which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"} i.e., the circumcenter of each triangle.
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the \code{Xp} points, Delaunay triangles based on \code{Yp} and also the AS proximity regions
#' for \code{Xp} points inside the convex hull of \code{Yp} points
#'
#' @seealso \code{\link{plotASregs.tri}}, \code{\link{plotPEregs.tri}}, \code{\link{plotPEregs}},
#' \code{\link{plotCSregs.tri}}, and \code{\link{plotCSregs}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-10 ; ny<-5
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3) #or M="CC"
#'
#' plotASregs(Xp,Yp,M,xlab="",ylab="")
#'
#' plotASregs(Xp,Yp[1:3,],M,xlab="",ylab="")
#'
#' Xp<-c(.5,.5)
#' plotASregs(Xp,Yp,M,xlab="",ylab="")
#' }
#'
#' @export plotASregs
plotASregs <- function(Xp,Yp,M="CC",main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
{
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

  if (nrow(Yp)==3)
  {
    plotASregs.tri(Xp,Yp,M,main,xlab,ylab,xlim,ylim)
  } else
  {
    if ((!is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates or "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)  #the Xp points inside the convex hull of Yp points

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)  #number of Delaunay triangles
    nx2<-nrow(Xch)  #number of Xp points inside the convex hull of Yp points

    if (is.null(xlim))
    {xlim<-range(Yp[,1],Xp[,1])
    xd<-xlim[2]-xlim[1]
    xlim<-xlim+xd*c(-.05,.05)}

    if (is.null(ylim))
    {ylim<-range(Yp[,2],Xp[,2])
    yd<-ylim[2]-ylim[1]
    ylim<-ylim+yd*c(-.05,.05)}

    if (is.null(main))
    {if (identical(M,"CC")) {
      main="AS Proximity Regions with CC-Vertex Regions"
    } else {Mvec= paste(M, collapse=",")
    main=paste("AS Proximity Regions with M = (",Mvec,")")}
    }

    plot(rbind(Xp),asp=1,main=main, xlab=xlab, ylab=ylab,
         xlim=xlim,ylim=ylim,pch=".",cex=3)#,...)

    if (nx2==0)
    {
      for (i in 1:nt)
      {
        tri<-Yp[DTr[i,],]  #vertices of the ith triangle
        polygon(tri,lty=2)
      }
    } else
    {
      i.tr<-rep(0,nx2)  #the vector of indices for the triangles that contain the Xp points
      for (i1 in 1:nx2)
        for (j1 in 1:nt)
        {
          Tri<-Yp[DTr[j1,],]
          if (in.triangle(Xch[i1,],Tri,boundary=TRUE)$in.tri )
            i.tr[i1]<-j1
        }


      for (i in 1:nt)
      {
        Tri<-Yp[DTr[i,],] #vertices of the ith triangle
        tri<-as.bas.tri(Tri)$tri #convert the triangle Tri into an nonscaled basic triangle, see as.bas.tri help page

        polygon(tri,lty=2)
        Xtri<-matrix(Xch[i.tr==i,],ncol=2)  #Xp points inside triangle i
        ni<-nrow(Xtri)
        if (ni>=1)
        {
          ################
          for (j in 1:ni)
          {
            P1<-Xtri[j,]
            rv<-ifelse(identical(M,"CC"), rv.triCC(P1,tri)$rv, rv.tri.cent(P1,tri,M)$rv)  #vertex region for P1
            RV<-tri[rv,]
            rad<-Dist(P1,RV)

            Int.Pts<-NAStri(P1,tri,M)
            L<-matrix(Int.Pts$L,ncol=2); R<-matrix(Int.Pts$R,ncol=2)
            segments(L[,1], L[,2], R[,1], R[,2], lty=1,col=2)
            Arcs<-Int.Pts$a;
            if (!is.null(Arcs))
            {
              K<-nrow(Arcs)/2
              for (k in 1:K)
              {A1<-Arcs[2*k-1,]; A2<-Arcs[2*k,];
              angles<-angle.str2end(A1,P1,A2)$c

              test.ang1<-angles[1]+(.01)*(angles[2]-angles[1])
              test.Pnt<-P1+rad*c(cos(test.ang1),sin(test.ang1))
              if (!in.triangle(test.Pnt,tri,boundary = TRUE)$i) {angles<-c(min(angles),max(angles)-2*pi)}
              plotrix::draw.arc(P1[1],P1[2],rad,angle1=angles[1],angle2=angles[2],col=2)
              }
            }
          }
          ################
        }
      }
    }
  }
} #end of the function
#'
