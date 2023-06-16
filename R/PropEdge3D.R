#PropEdge3D.R;

#################################################################
#Functions for NPE in R^3
#################################################################

#' @title The vertices of the Proportional Edge (PE) Proximity Region in the standard regular tetrahedron
#'
#' @description Returns the vertices of the PE proximity region (which is itself a tetrahedron) for a point in the
#' standard regular tetrahedron \eqn{T_h=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))=}
#' \code{(rv=1,rv=2,rv=3,rv=4)}.
#'
#' PE proximity region is defined with respect to the tetrahedron \eqn{T_h}
#' with expansion parameter \eqn{r \ge 1} and vertex regions based on the circumcenter of \eqn{T_h} (which is equivalent
#' to the center of mass in the standard regular tetrahedron).
#'
#' Vertex regions are labeled as \code{1,2,3,4} rowwise for the vertices of the tetrahedron \eqn{T_h}.
#' \code{rv} is the index of the vertex region \code{p} resides, with default=\code{NULL}.
#' If \code{p} is outside of \eqn{T_h}, it returns \code{NULL} for the proximity region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p A 3D point whose PE proximity region is to be computed.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param rv Index of the vertex region containing the point, either \code{1,2,3,4} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return Vertices of the tetrahedron which constitutes the PE proximity region with expansion parameter
#' \code{r} and circumcenter (or center of mass) for a point \code{p} in the standard regular tetrahedron
#'
#' @seealso \code{\link{NPEtetra}}, \code{\link{NPEtri}} and \code{\link{NPEint}}
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
#' n<-3
#' Xp<-runif.std.tetra(n)$g
#' r<-1.5
#' NPEstd.tetra(Xp[1,],r)
#'
#' #or try
#' RV<-rel.vert.tetraCC(Xp[1,],tetra)$rv
#' NPEstd.tetra(Xp[1,],r,rv=RV)
#'
#' NPEstd.tetra(c(-1,-1,-1),r,rv=NULL)
#' }
#'
#' @export NPEstd.tetra
NPEstd.tetra <- function(p,r,rv=NULL)
{
  if (!is.point(p,3) )
  {stop('p must be a numeric 3D point')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  th<-rbind(A,B,C,D)  #standard regular tetrahedron

  if (!in.tetrahedron(p,th,boundary=TRUE)$in.tetra)
  {reg<-NULL; return(reg); stop}

  if (is.null(rv))
  {rv<-rel.vert.tetraCC(p,th)$rv #vertex region for p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3,4))!=1)
  {stop('vertex index, rv, must be 1, 2, 3 or 4')}}

  p1<-p[1]; p2<-p[2]; p3<-p[3];

  if (rv==1)
  {
    A1<-c((1/6*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2+6*p1))*r, 0, 0);
    A2<-c((1/12*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2+6*p1))*r, (1/12)*sqrt(3)*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2+6*p1)*r, 0);
    A3<-c((1/12*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2+6*p1))*r, (1/36)*sqrt(3)*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2+6*p1)*r, (1/18)*sqrt(3)*sqrt(2)*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2+6*p1)*r);
    reg<-rbind(A,A1,A2,A3)
  } else if (rv==2)
  {
    B1<-c((1/2)*p1*r-(1/2)*r-(1/6)*sqrt(3)*p2*r+1-(1/12)*sqrt(3)*sqrt(2)*p3*r, (1/12)*sqrt(3)*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2-6*p1+6)*r, 0);
    B2<-c((1/2)*p1*r-(1/2)*r-(1/6)*sqrt(3)*p2*r+1-(1/12)*sqrt(3)*sqrt(2)*p3*r, (1/36)*sqrt(3)*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2-6*p1+6)*r, (1/18)*sqrt(3)*sqrt(2)*(sqrt(3)*sqrt(2)*p3+2*sqrt(3)*p2-6*p1+6)*r);
    B3<-c(p1*r-r-(1/3)*sqrt(3)*p2*r+1-(1/6)*sqrt(3)*sqrt(2)*p3*r, 0, 0);
    reg<-rbind(B,B1,B2,B3)
  } else if (rv==3)
  {
    C1<-c(-(1/12)*sqrt(3)*(sqrt(2)*p3*r+2*sqrt(3)*r-4*p2*r-2*sqrt(3)), -(1/4)*sqrt(2)*p3*r-(1/2)*sqrt(3)*r+p2*r+(1/2)*sqrt(3), 0);
    C2<-c(1/2, -(1/3)*sqrt(3)*r+(2/3)*p2*r+(1/2)*sqrt(3)-(1/6)*sqrt(2)*p3*r, (1/3)*sqrt(3)*sqrt(2)*r-(2/3)*sqrt(2)*p2*r+(1/3)*p3*r);
    C3<-c(-(1/12)*sqrt(3)*(-sqrt(2)*p3*r-2*sqrt(3)*r+4*p2*r-2*sqrt(3)), -(1/4)*sqrt(2)*(p3*r+sqrt(3)*sqrt(2)*r-2*sqrt(2)*p2*r-sqrt(3)*sqrt(2)), 0);
    reg<-rbind(C,C1,C2,C3)
  } else {
    D1<-c(-(1/12)*sqrt(6)*(sqrt(6)*r-3*p3*r-sqrt(6)), -(1/12)*sqrt(2)*(sqrt(3)*sqrt(2)*r-sqrt(3)*sqrt(2)-3*p3*r), -(1/3)*sqrt(3)*sqrt(2)*r+(1/3)*sqrt(3)*sqrt(2)+p3*r);
    D2<-c((1/12*(sqrt(6)*r-3*p3*r+sqrt(6)))*sqrt(6), -(1/12)*sqrt(2)*(sqrt(3)*sqrt(2)*r-sqrt(3)*sqrt(2)-3*p3*r), -(1/3)*sqrt(3)*sqrt(2)*r+(1/3)*sqrt(3)*sqrt(2)+p3*r);
    D3<-c(1/2, -(1/36)*sqrt(6)*sqrt(3)*(6*p3*r-2*r*sqrt(6)-sqrt(6)), -(1/3)*r*sqrt(6)+p3*r+(1/3)*sqrt(6));
    reg<-rbind(D,D1,D2,D3)
  }

  vec1<-rep(1,4);
  D0<-abs(det(matrix(cbind(th,vec1),ncol=4)))
  D1<-abs(det(matrix(cbind(reg,vec1),ncol=4)))
  if (D1>=D0)
  {reg<-th}
  row.names(reg)<-c()
  reg
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Proportional Edge Proximity Catch
#' Digraphs (PE-PCDs) - standard regular tetrahedron case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))} for points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' returns 0 otherwise, where \eqn{N_{PE}(x,r)} is the PE proximity region for point \eqn{x} with expansion parameter \eqn{r \ge 1}.
#'
#' PE proximity region is defined with respect to the standard regular tetrahedron
#' \eqn{T_h=T(v=1,v=2,v=3,v=4)=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))} and vertex regions
#' are based on the circumcenter (which is equivalent to the center of mass for standard regular tetrahedron)
#' of \eqn{T_h}. \code{rv} is the index of the vertex region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_h}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1 A 3D point whose PE proximity region is constructed.
#' @param p2 A 3D point. The function determines whether \code{p2} is inside the PE proximity region of
#' \code{p1} or not.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param rv Index of the vertex region containing the point, either \code{1,2,3,4} (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))} for points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' returns 0 otherwise
#'
#' @seealso \code{\link{IarcPEtetra}}, \code{\link{IarcPEtri}} and \code{\link{IarcPEint}}
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
#' n<-3  #try also n<-20
#' Xp<-runif.std.tetra(n)$g
#' r<-1.5
#' IarcPEstd.tetra(Xp[1,],Xp[3,],r)
#' IarcPEstd.tetra(c(.4,.4,.4),c(.5,.5,.5),r)
#'
#' #or try
#' RV<-rel.vert.tetraCC(Xp[1,],tetra)$rv
#' IarcPEstd.tetra(Xp[1,],Xp[3,],r,rv=RV)
#'
#' P1<-c(.1,.1,.1)
#' P2<-c(.5,.5,.5)
#' IarcPEstd.tetra(P1,P2,r)
#' }
#'
#' @export IarcPEstd.tetra
IarcPEstd.tetra <- function(p1,p2,r,rv=NULL)
{
  if (!is.point(p1,3) || !is.point(p2,3))
  {stop('p1 and p2 must be a numeric 3D point')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  th<-rbind(A,B,C,D)  #standard regular tetrahedron

  if (!in.tetrahedron(p1,th,boundary=TRUE)$in.tetra || !in.tetrahedron(p2,th,boundary=TRUE)$in.tetra)
  {arc<-0
  } else
  {
    if (is.null(rv))
    {rv<-rel.vert.tetraCC(p1,th)$rv #vertex region for p1
    } else
    {  if (!is.numeric(rv) || sum(rv==c(1,2,3,4))!=1)
    {stop('vertex index, rv, must be 1, 2, 3 or 4')}}

    y1<-th[1,]; y2<-th[2,]; y3<-th[3,]; y4<-th[4,];
    X1<-p1[1]; Y1<-p1[2]; Z1<-p1[3];
    X2<-p2[1]; Y2<-p2[2]; Z2<-p2[3];
    arc<-0;
    if (rv==1)
    {
      x1n<-X1*r; y1n<-Y1*r; z1n<-Z1*r;
      if ( Z2 < as.vector(paraplane(c(x1n,y1n,z1n),y2,y3,y4,X2,Y2)$z) ) {arc <-1}
    } else {
      if (rv==2)
      {
        x1n<-1+(X1-1)*r; y1n<-Y1*r;  z1n<-Z1*r;
        if ( Z2 < as.vector(paraplane(c(x1n,y1n,z1n),y1,y3,y4,X2,Y2)$z) ) {arc <-1}
      } else {
        if (rv==3)
        {
          x1n<-1/2+(X1-1/2)*r; y1n<-sqrt(3)/2+(Y1-sqrt(3)/2)*r; z1n<-Z1*r;
          if ( Z2 < as.vector(paraplane(c(x1n,y1n,z1n),y1,y2,y4,X2,Y2)$z) ) {arc<-1}
        } else {
          z1n<-y4[3]+(Z1-y4[3])*r;
          if ( Z2 > z1n ) {arc<-1}
        }}}
  }
  arc
} #end of the function
#'

#################################################################

#' @title The plot of the Proportional Edge (PE) Proximity Regions for a 3D data set - standard
#' regular tetrahedron case
#'
#' @description Plots the points in and outside of the standard regular tetrahedron
#' \eqn{T_h=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))} and also the PE proximity regions
#' for points in data set \code{Xp}.
#'
#' PE proximity regions are defined with respect to the standard regular tetrahedron \eqn{T_h}
#' with expansion parameter \eqn{r \ge 1}, so PE proximity regions are defined only for points inside \eqn{T_h}.
#'
#' Vertex regions are based on circumcenter (which is equivalent to the center of mass for the standard
#' regular tetrahedron) of \eqn{T_h}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 3D points for which PE proximity regions are constructed.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab,zlab titles for the \eqn{x}, \eqn{y}, and \eqn{z} axes, respectively (default=\code{NULL} for all).
#' @param xlim,ylim,zlim Two \code{numeric} vectors of length 2, giving the \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinate ranges
#' (default=\code{NULL} for all).
#' @param \dots Additional \code{scatter3D} parameters.
#'
#' @return Plot of the PE proximity regions for points inside the standard regular tetrahedron \eqn{T_h}
#' (and just the points outside \eqn{T_h})
#'
#' @seealso \code{\link{plotPEregs}}, \code{\link{plotASregs.tri}}, \code{\link{plotASregs}},
#' \code{\link{plotCSregs.tri}}, and \code{\link{plotCSregs}}
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
#' r<-1.5
#'
#' n<-3  #try also n<-20
#' Xp<-runif.std.tetra(n)$g  #try also Xp[,1]<-Xp[,1]+1
#'
#' plotPEregs.std.tetra(Xp[1:3,],r)
#'
#' P1<-c(.1,.1,.1)
#' plotPEregs.std.tetra(rbind(P1,P1),r)
#' }
#'
#' @export plotPEregs.std.tetra
plotPEregs.std.tetra <- function(Xp,r,main=NULL,xlab=NULL,ylab=NULL,zlab=NULL,xlim=NULL,ylim=NULL,zlim=NULL, ...)
{
  if (!is.numeric(as.matrix(Xp)) )
  {stop('Xp must be numeric')}

  if (is.point(Xp,dim=length(Xp)))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  n<-nrow(Xp)

  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  th<-rbind(A,B,C,D)  #standard regular tetrahedron

  in.tetra<-rep(0,n)
  for (i in 1:n)
    in.tetra[i]<-in.tetrahedron(Xp[i,],th,boundary=TRUE)$in.tetra #indices of the Xp points inside the tetrahedron

  Xtetra<-matrix(Xp[in.tetra==1,],ncol=3)  #the Xp points inside the tetrahedron
  nt<-length(Xtetra)/3 #number of Xp points inside the tetrahedron

  if (is.null(xlim))
  {xlim<-range(th[,1],Xp[,1])}
  if (is.null(ylim))
  {ylim<-range(th[,2],Xp[,2])}
  if (is.null(zlim))
  {zlim<-range(th[,3],Xp[,3])}

  xr<-xlim[2]-xlim[1]
  yr<-ylim[2]-ylim[1]
  zr<-zlim[2]-zlim[1]

  if (is.null(main))
  { main=paste("PE Proximity Regions with r = ",r,"\n for Data in Standard Regular Tetrahedron",sep="")}

  plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi =0,theta=40, bty = "g",main=main,xlab=xlab, ylab=ylab, zlab=zlab,
                    xlim=xlim+xr*c(-.05,.05), ylim=ylim+yr*c(-.05,.05), zlim=zlim+zr*c(-.05,.05),
                    pch = 20, cex = 1, ticktype = "detailed", ...)
  #add the vertices of the tetrahedron
  plot3D::points3D(th[,1],th[,2],th[,3], add = TRUE)
  L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
  plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=1,lty=2)
  if (nt>=1)
  {
    for (i in 1:nt)
    {
      P1<-as.vector(Xtetra[i,])
      RV<-rel.vert.tetraCC(P1,th)$rv

      pr<-NPEstd.tetra(P1,r,rv=RV)
      L<-rbind(pr[1,],pr[1,],pr[1,],pr[2,],pr[2,],pr[3,]);
      R<-rbind(pr[2,],pr[3,],pr[4,],pr[3,],pr[4,],pr[4,])
      plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2,col="blue")
    }
  }
} #end of the function
#'

#################################################################

#' @title The vertices of the Proportional Edge (PE) Proximity Region in a tetrahedron
#'
#' @description Returns the vertices of the PE proximity region (which is itself a tetrahedron) for a point in the
#' tetrahedron \code{th}.
#'
#' PE proximity region is defined with respect to the tetrahedron \code{th}
#' with expansion parameter \eqn{r \ge 1} and vertex regions based on the center \code{M} which is circumcenter (\code{"CC"}) or
#' center of mass (\code{"CM"}) of \code{th} with default=\code{"CM"}.
#'
#' Vertex regions are labeled as \code{1,2,3,4} rowwise for the vertices of the tetrahedron \code{th}.
#' \code{rv} is the index of the vertex region \code{p} resides, with default=\code{NULL}.
#' If \code{p} is outside of \code{th}, it returns \code{NULL} for the proximity region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p A 3D point whose PE proximity region is to be computed.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param rv Index of the vertex region containing the point, either \code{1,2,3,4} (default is \code{NULL}).
#'
#' @return Vertices of the tetrahedron which constitutes the PE proximity region with expansion parameter
#' \code{r} and circumcenter (or center of mass) for a point \code{p} in the tetrahedron
#'
#' @seealso \code{\link{NPEstd.tetra}}, \code{\link{NPEtri}} and \code{\link{NPEint}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' set.seed(1)
#' tetra<-rbind(A,B,C,D)+matrix(runif(12,-.25,.25),ncol=3)
#' n<-3  #try also n<-20
#'
#' Xp<-runif.tetra(n,tetra)$g
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.5
#'
#' NPEtetra(Xp[1,],tetra,r)  #uses the default M="CM"
#' NPEtetra(Xp[1,],tetra,r,M="CC")
#'
#' #or try
#' RV<-rel.vert.tetraCM(Xp[1,],tetra)$rv
#' NPEtetra(Xp[1,],tetra,r,M,rv=RV)
#'
#' P1<-c(.1,.1,.1)
#' NPEtetra(P1,tetra,r,M)
#' }
#'
#' @export NPEtetra
NPEtetra <- function(p,th,r,M="CM",rv=NULL)
{
  if (!is.point(p,3) )
  {stop('p must be a numeric 3D point')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  if (!in.tetrahedron(p,th,boundary=T)$in.tetra)
  {reg<-NULL; return(reg); stop}

  Rv<-rv
  if (is.null(Rv))
  {Rv<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p,th)$rv,rel.vert.tetraCM(p,th)$rv)  #vertex region for pt
  } else
  {  if (!is.numeric(Rv) || sum(Rv==c(1,2,3,4))!=1)
  {stop('vertex index, rv, must be 1, 2, 3 or 4')}}

  A<-th[Rv,];
  th0<-th[-Rv,]; B<-th0[1,]; C<-th0[2,]; D<-th0[3,]
  d1<-dist.point2plane(A,B,C,D)$d; d2<-dist.point2plane(p,B,C,D)$d;
  rd<-(d1-d2)/d1; #distance ratio
  A1p<-A+rd*(B-A); A1<-A+r*(A1p-A);
  A2p<-A+rd*(C-A); A2<-A+r*(A2p-A);
  A3p<-A+rd*(D-A); A3<-A+r*(A3p-A);
  reg<-rbind(A,A1,A2,A3)

  vec1<-rep(1,4);
  D0<-abs(det(matrix(cbind(th,vec1),ncol=4)))
  D1<-abs(det(matrix(cbind(reg,vec1),ncol=4)))
  if (D1>=D0)
  {reg<-th}
  row.names(reg)<-c()
  reg
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from one 3D point to another 3D point for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))} for 3D points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' returns 0 otherwise, where N_{PE}(x,r) is the PE proximity region for point \eqn{x} with the expansion parameter \eqn{r \ge 1}.
#'
#' PE proximity region is constructed with respect to the tetrahedron \code{th} and
#' vertex regions are based on the center \code{M} which is circumcenter (\code{"CC"}) or
#' center of mass (\code{"CM"}) of \code{th} with default=\code{"CM"}.
#' \code{rv} is the index of the vertex region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \code{th}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1 A 3D point whose PE proximity region is constructed.
#' @param p2 A 3D point. The function determines whether \code{p2} is inside the PE proximity region of
#' \code{p1} or not.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param rv Index of the \code{M}-vertex region containing the point, either \code{1,2,3,4}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))} for \code{p1}, that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)}, returns 0 otherwise
#'
#' @seealso \code{\link{IarcPEstd.tetra}}, \code{\link{IarcPEtri}} and \code{\link{IarcPEint}}
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
#' n<-3  #try also n<-20
#'
#' Xp<-runif.tetra(n,tetra)$g
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.5
#'
#' IarcPEtetra(Xp[1,],Xp[2,],tetra,r)  #uses the default M="CM"
#' IarcPEtetra(Xp[1,],Xp[2,],tetra,r,M)
#'
#' IarcPEtetra(c(.4,.4,.4),c(.5,.5,.5),tetra,r,M)
#'
#' #or try
#' RV<-rel.vert.tetraCC(Xp[1,],tetra)$rv
#' IarcPEtetra(Xp[1,],Xp[3,],tetra,r,M,rv=RV)
#'
#' P1<-c(.1,.1,.1)
#' P2<-c(.5,.5,.5)
#' IarcPEtetra(P1,P2,tetra,r,M)
#' }
#'
#' @export IarcPEtetra
IarcPEtetra <- function(p1,p2,th,r,M="CM",rv=NULL)
{
  if (!is.point(p1,3) || !is.point(p2,3))
  {stop('p1 and p2 must be a numeric 3D point')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  if (!in.tetrahedron(p1,th,boundary=TRUE)$in.tetra || !in.tetrahedron(p2,th,boundary=TRUE)$in.tetra)
  {arc<-0
  } else
  {
    if (is.null(rv))
    {rv<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p1,th)$rv,rel.vert.tetraCM(p1,th)$rv)  #vertex region for p1
    } else
    {  if (!is.numeric(rv) || sum(rv==c(1,2,3,4))!=1)
    {stop('vertex index, rv, must be 1, 2, 3 or 4')}}

    pr<-NPEtetra(p1,th,r,M,rv)
    D0<-det(matrix(cbind(pr,rep(1,4)),ncol=4))
    arc<-ifelse(isTRUE(all.equal(D0, 0)),0,sum(in.tetrahedron(p2,pr,boundary = TRUE)$i))
  }
  arc
} #end of the function
#'

################################################################

#' @title Incidence matrix for Proportional Edge Proximity Catch Digraphs (PE-PCDs) - one tetrahedron case
#'
#' @description Returns the incidence matrix for the PE-PCD whose vertices are the given 3D numerical data set, \code{Xp},
#' in the tetrahedron \eqn{th=T(v=1,v=2,v=3,v=4)}.
#'
#' PE proximity regions are constructed with respect to tetrahedron
#' \code{th} with expansion parameter \eqn{r \ge 1} and vertex regions are based on the center \code{M} which is circumcenter (\code{"CC"})
#' or center of mass (\code{"CM"}) of \code{th} with default=\code{"CM"}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 3D points which constitute the vertices of PE-PCD.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#'
#' @return Incidence matrix for the PE-PCD with vertices being 3D data set, \code{Xp},
#' in the tetrahedron \code{th} with vertex regions based on circumcenter or center of mass
#'
#' @seealso \code{\link{inci.matPEtri}}, \code{\link{inci.matPE1D}}, and  \code{\link{inci.matPE}}
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
#' n<-5
#'
#' Xp<-runif.tetra(n,tetra)$g  #try also Xp<-c(.5,.5,.5)
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.5
#'
#' IM<-inci.matPEtetra(Xp,tetra,r=1.25)  #uses the default M="CM"
#' IM<-inci.matPEtetra(Xp,tetra,r=1.25,M)
#' IM
#' dom.num.greedy(IM)
#' Idom.num.up.bnd(IM,3) #try also dom.num.exact(IM)  #this might take a long time for large n
#' }
#'
#' @export inci.matPEtetra
inci.matPEtetra <- function(Xp,th,r,M="CM")
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp,dim=length(Xp)))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  n<-nrow(Xp)

  inc.mat<-matrix(0, nrow=n, ncol=n)
  if (n>1)
  {
    for (i in 1:n)
    {p1<-Xp[i,]
    v<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p1,th)$rv,rel.vert.tetraCM(p1,th)$rv)  #vertex region for p1
    for (j in ((1:n)) )
    {p2<-Xp[j,]
    inc.mat[i,j]<-IarcPEtetra(p1,p2,th,r,M,rv=v)
    }
    }
  }
  diag(inc.mat)<-1
  inc.mat
} #end of the function
#'

#################################################################

#' @title Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' and quantities related to the tetrahedron - one tetrahedron case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' whose vertices are the
#' given 3D numerical data set, \code{Xp}.
#' It also provides number of vertices
#' (i.e., number of data points inside the tetrahedron)
#' and indices of the data points that reside in the tetrahedron.
#'
#' PE proximity region is constructed with respect to the tetrahedron \code{th} and
#' vertex regions are based on the center \code{M} which is circumcenter (\code{"CC"}) or
#' center of mass (\code{"CM"}) of \code{th} with default=\code{"CM"}.
#' For the number of arcs, loops are not
#' allowed so arcs are only possible for points inside the tetrahedron \code{th} for this function.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 3D points which constitute the vertices of PE-PCD.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and quantities related to the tetrahedron}
#' \item{num.arcs}{Number of arcs of the PE-PCD}
#' \item{num.in.tetra}{Number of \code{Xp} points in the tetrahedron, \code{th}}
#' \item{ind.in.tetra}{The vector of indices of the \code{Xp} points that reside in the tetrahedron}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support tetrahedron.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsPEtri}}, \code{\link{num.arcsCStri}},
#' and \code{\link{num.arcsAStri}}
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
#' n<-10  #try also n<-20
#' set.seed(1)
#' Xp<-runif.tetra(n,tetra)$g
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.25
#'
#' Narcs = num.arcsPEtetra(Xp,tetra,r,M)
#' Narcs
#' summary(Narcs)
#' #plot(Narcs)
#' }
#'
#' @export num.arcsPEtetra
num.arcsPEtetra <- function(Xp,th,r,M="CM")
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('The tetrahedron th is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  n<-nrow(Xp)
  arcs<-0
  ind.in.tetra = NULL
  if (n<=0)
  {
    arcs<-0
  } else
  {
    for (i in 1:n)
    { p1 = Xp[i,]
    if (in.tetrahedron(p1,th,boundary=TRUE)$in.tetra)
    {  vert<- ifelse(identical(M,"CC"),
                     rel.vert.tetraCC(p1,th)$rv,
                     rel.vert.tetraCM(p1,th)$rv)  #vertex region for p1
    ind.in.tetra = c(ind.in.tetra,i)
    for (j in (1:n)[-i])  #to avoid loops
    {  p2 = Xp[j,]
    arcs<-arcs+IarcPEtetra(p1,p2,th,r,M,rv=vert)
    }
    }
    }
  }

  NinTetra = length(ind.in.tetra)
  desc<-"Number of Arcs of the PE-PCD with vertices Xp and Quantities Related to the Support Tetrahedron"

  res<-list(desc=desc, #description of the output
            num.arcs=arcs, #number of arcs for the PE-PCD
            num.in.tetra=NinTetra, # number of Xp points in CH of Yp points
            ind.in.tetra=ind.in.tetra, #indices of data points inside the tetrahedron
            tess.points=th, #tessellation points
            vertices=Xp #vertices of the digraph
  )

  class(res)<-"NumArcs"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Arc density of Proportional Edge Proximity Catch Digraphs (PE-PCDs) - one tetrahedron case
#'
#' @description Returns the arc density of PE-PCD whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the tetrahedron \code{th}.
#'
#' PE proximity region is constructed with respect to the tetrahedron \code{th} and
#' vertex regions are based on the center \code{M} which is circumcenter (\code{"CC"}) or
#' center of mass (\code{"CM"}) of \code{th} with default=\code{"CM"}.
#' For the number of arcs, loops are not
#' allowed so arcs are only possible for points inside the tetrahedron \code{th} for this function.
#'
#' \code{th.cor} is a logical argument for tetrahedron correction (default is \code{TRUE}), if \code{TRUE}, only the points
#' inside the tetrahedron are considered (i.e., digraph induced by these vertices are considered) in computing
#' the arc density, otherwise all points are considered (for the number of vertices in the denominator of arc
#' density).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the PE-PCD.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param th.cor A logical argument for computing the arc density for only the points inside the tetrahedron,
#' \code{th}. (default is \code{th.cor=FALSE}), i.e., if \code{th.cor=TRUE} only the induced digraph
#' with the vertices inside \code{th} are considered in the computation of arc density.
#'
#' @return Arc density of PE-PCD whose vertices are the 2D numerical data set, \code{Xp};
#' PE proximity regions are defined with respect to the tetrahedron \code{th} and \code{M}-vertex regions
#'
#' @seealso \code{\link{PEarc.dens.tri}} and \code{\link{num.arcsPEtetra}}
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
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.tetra(n,tetra)$g
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.5
#'
#' num.arcsPEtetra(Xp,tetra,r,M)
#' PEarc.dens.tetra(Xp,tetra,r,M)
#' PEarc.dens.tetra(Xp,tetra,r,M,th.cor = FALSE)
#' }
#'
#' @export PEarc.dens.tetra
PEarc.dens.tetra <- function(Xp,th,r,M="CM",th.cor=FALSE)
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron th is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  nx<-nrow(Xp)
  narcs<-num.arcsPEtetra(Xp,th,r,M)$num.arcs
  # mean.rho<-muPE2D(r)
  # var.rho<-asyvarPE2D(r)

  if (th.cor==TRUE)
  {
    ind.it<-c()
    for (i in 1:nx)
    {
      ind.it<-c(ind.it,in.tetrahedron(Xp[i,],th)$in.tetra)
    }
    dat.it<-Xp[ind.it,] #Xp points inside the tetrahedron
    NinTetra<-nrow(dat.it)
    if (NinTetra<=1)
    {stop('There are not enough Xp points in the tetrahedron, th, to compute the arc density!')}
    n<-NinTetra
  } else
  {
    n<-nx
  }
  rho<-narcs/(n*(n-1))
  # std.rho<-sqrt(n)*(rho-mean.rho)/sqrt(var.rho)
  # list(
  rho#arc.dens=rho #, #arc density
  #$   std.arc.dens=std.rho, #standardized arc density
  #   caveat="The standardized arc density is only correct when M is the center of mass in the current version." #caveat
  # )
} #end of the function
#'

#################################################################

#' @title The plot of the Proportional Edge (PE) Proximity Regions for a 3D data set - one tetrahedron case
#'
#' @description Plots the points in and outside of the tetrahedron \code{th} and also the PE proximity regions  (which are also
#' tetrahedrons) for points inside the tetrahedron \code{th}.
#'
#' PE proximity regions are constructed with respect to
#' tetrahedron \code{th} with expansion parameter \eqn{r \ge 1} and vertex regions are based on the center \code{M} which is
#' circumcenter (\code{"CC"}) or center of mass (\code{"CM"}) of \code{th} with default=\code{"CM"}, so PE proximity regions are defined
#' only for points inside the tetrahedron \code{th}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 3D points for which PE proximity regions are constructed.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab,zlab Titles for the \eqn{x}, \eqn{y}, and \eqn{z} axes, respectively (default=\code{NULL} for all).
#' @param xlim,ylim,zlim Two \code{numeric} vectors of length 2, giving the \eqn{x}-, \eqn{y}-, and \eqn{z}-coordinate ranges
#' (default=\code{NULL} for all).
#' @param \dots Additional \code{scatter3D} parameters.
#'
#' @return Plot of the PE proximity regions for points inside the tetrahedron \code{th}
#' (and just the points outside \code{th})
#'
#' @seealso \code{\link{plotPEregs.std.tetra}}, \code{\link{plotPEregs.tri}} and \code{\link{plotPEregs.int}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' set.seed(1)
#' tetra<-rbind(A,B,C,D)+matrix(runif(12,-.25,.25),ncol=3) #adding jitter to make it non-regular
#'
#' n<-5  #try also n<-20
#' Xp<-runif.tetra(n,tetra)$g  #try also Xp[,1]<-Xp[,1]+1
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.5
#'
#' plotPEregs.tetra(Xp,tetra,r) #uses the default M="CM"
#' plotPEregs.tetra(Xp,tetra,r,M="CC")
#'
#' plotPEregs.tetra(Xp[1,],tetra,r)  #uses the default M="CM"
#' plotPEregs.tetra(Xp[1,],tetra,r,M)
#' }
#'
#' @export plotPEregs.tetra
plotPEregs.tetra <- function(Xp,th,r,M="CM",main=NULL,xlab=NULL,ylab=NULL,zlab=NULL,xlim=NULL,ylim=NULL,zlim=NULL, ...)
{
  if (!is.numeric(as.matrix(Xp)) )
  {stop('Xp must be numeric')}

  if (is.point(Xp,dim=length(Xp)))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  n<-nrow(Xp)

  in.tetra<-rep(0,n)
  for (i in 1:n)
    in.tetra[i]<-in.tetrahedron(Xp[i,],th,boundary=TRUE)$in.tetra #indices of the Xp points inside the tetrahedron

  Xtetra<-matrix(Xp[in.tetra==1,],ncol=3)  #the Xp points inside the tetrahedron
  nt<-length(Xtetra)/3 #number of Xp points inside the tetrahedron

  if (is.null(xlim))
  {xlim<-range(th[,1],Xp[,1])}
  if (is.null(ylim))
  {ylim<-range(th[,2],Xp[,2])}
  if (is.null(zlim))
  {zlim<-range(th[,3],Xp[,3])}

  xr<-xlim[2]-xlim[1]
  yr<-ylim[2]-ylim[1]
  zr<-zlim[2]-zlim[1]

  if (is.null(main))
  { main=paste("PE Proximity Regions with r = ",r," and M = ",M,sep="")}

  plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi =0,theta=40, bty = "g",main=main,xlab=xlab, ylab=ylab, zlab=zlab,
                    xlim=xlim+xr*c(-.05,.05), ylim=ylim+yr*c(-.05,.05), zlim=zlim+zr*c(-.05,.05),
                    pch = 20, cex = 1, ticktype = "detailed", ...)
  #add the vertices of the tetrahedron
  plot3D::points3D(th[,1],th[,2],th[,3], add = TRUE)
  A<-th[1,]; B<-th[2,]; C<-th[3,]; D<-th[4,]
  L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
  plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=1,lty=2)
  if (nt>=1)
  {
    for (i in 1:nt)
    {
      P1<-as.vector(Xtetra[i,])
      RV<-ifelse(identical(M,"CC"),rel.vert.tetraCC(P1,th)$rv,rel.vert.tetraCM(P1,th)$rv)
      pr<-NPEtetra(P1,th,r,M,rv=RV)
      L<-rbind(pr[1,],pr[1,],pr[1,],pr[2,],pr[2,],pr[3,]);
      R<-rbind(pr[2,],pr[3,],pr[4,],pr[3,],pr[4,],pr[4,])
      plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2,col="blue")
    }
  }
} #end of the function
#'

#################################################################

#' @title The indicator for a 3D point being a dominating point for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - standard regular tetrahedron case
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp} in the
#' standard regular tetrahedron \eqn{T_h=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))}, that is,
#' returns 1 if \code{p} is a dominating point of PE-PCD, returns 0 otherwise.
#'
#' Point, \code{p}, is in the vertex region of vertex \code{rv} (default is \code{NULL}); vertices are labeled as \code{1,2,3,4}
#' in the order they are stacked row-wise in \eqn{T_h}.
#'
#' PE proximity region is constructed with respect to the tetrahedron \eqn{T_h} with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center of mass \eqn{CM} (equivalent to circumcenter in this case).
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p A 3D point that is to be tested for being a dominating point or not of the PE-PCD.
#' @param Xp A set of 3D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param rv Index of the vertex whose region contains point \code{p}, \code{rv} takes the vertex labels
#' as \code{1,2,3,4} as in the row order of the vertices in standard regular tetrahedron, default is \code{NULL}.
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num1PEtetra}}, \code{\link{Idom.num1PEtri}} and \code{\link{Idom.num1PEbasic.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' tetra<-rbind(A,B,C,D)
#'
#' n<-5 #try also n<-20
#' Xp<-runif.std.tetra(n)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))
#' r<-1.5
#'
#' P<-c(.4,.1,.2)
#' Idom.num1PEstd.tetra(Xp[1,],Xp,r)
#' Idom.num1PEstd.tetra(P,Xp,r)
#'
#' Idom.num1PEstd.tetra(Xp[1,],Xp,r)
#' Idom.num1PEstd.tetra(Xp[1,],Xp[1,],r)
#'
#' #or try
#' RV<-rel.vert.tetraCC(Xp[1,],tetra)$rv
#' Idom.num1PEstd.tetra(Xp[1,],Xp,r,rv=RV)
#'
#' Idom.num1PEstd.tetra(c(-1,-1,-1),Xp,r)
#' Idom.num1PEstd.tetra(c(-1,-1,-1),c(-1,-1,-1),r)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1PEstd.tetra(Xp[i,],Xp,r))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#' g1.pts<-Xp[ind.gam1,]
#'
#' Xlim<-range(tetra[,1],Xp[,1])
#' Ylim<-range(tetra[,2],Xp[,2])
#' Zlim<-range(tetra[,3],Xp[,3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi =0,theta=40, bty = "g",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
#'          pch = 20, cex = 1, ticktype = "detailed")
#' #add the vertices of the tetrahedron
#' plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#' if (length(g1.pts)!=0)
#' {
#'   if (length(g1.pts)==3) g1.pts<-matrix(g1.pts,nrow=1)
#'   plot3D::points3D(g1.pts[,1],g1.pts[,2],g1.pts[,3], pch=4,col="red", add=TRUE)}
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#'
#' CM<-apply(tetra,2,mean)
#' D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
#' L<-rbind(D1,D2,D3,D4,D5,D6); R<-matrix(rep(CM,6),ncol=3,byrow=TRUE)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
#'
#' P<-c(.4,.1,.2)
#' Idom.num1PEstd.tetra(P,Xp,r)
#'
#' Idom.num1PEstd.tetra(c(-1,-1,-1),Xp,r,ch.data.pnt = FALSE)
#' #gives an error message if ch.data.pnt = TRUE
#' }
#'
#' @export
Idom.num1PEstd.tetra <- function(p,Xp,r,rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p,3))
  {stop('p must be a numeric 3D point')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp,3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('point, p, is not a data point in Xp')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (isTRUE(all.equal(matrix(p,ncol=3),matrix(Xp,ncol=3))))
  {dom<-1;return(dom);stop}

  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  th<-rbind(A,B,C,D)  #standard regular tetrahedron
  if (in.tetrahedron(p,th)$in.tetra==FALSE)
  {dom<-0; return(dom); stop}

  if (is.null(rv))
  {rv<-rel.vert.tetraCC(p,th)$rv #vertex region for p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3,4))!=1)
  {stop('vertex index, rv, must be 1, 2, 3 or 4')}}

  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (IarcPEstd.tetra(p,Xp[i,],r,rv=rv)==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for two 3D points constituting a dominating set for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - standard regular tetrahedron case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp} in the standard regular
#' tetrahedron \eqn{T_h=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))}, that is,
#' returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD, returns 0 otherwise.
#'
#' Point, \code{p1}, is in the region of vertex \code{rv1} (default is \code{NULL}) and point, \code{p2}, is in the region of vertex \code{rv2}
#' (default is \code{NULL}); vertices (and hence \code{rv1} and \code{rv2}) are labeled as \code{1,2,3,4} in the order they are stacked
#' row-wise in \eqn{T_h}.
#'
#' PE proximity region is constructed with respect to the tetrahedron \eqn{T_h} with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center of mass \eqn{CM} (equivalent to circumcenter in this case).
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1} and \code{p2} are data points in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the points \code{p1} and \code{p2} would constitute a dominating set
#' if they actually were both in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1,p2 Two 3D points to be tested for constituting a dominating set of the PE-PCD.
#' @param Xp A set of 3D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param rv1,rv2 The indices of the vertices whose regions contains \code{p1} and \code{p2}, respectively.
#' They take the vertex labels as \code{1,2,3,4} as in the row order of the vertices in \eqn{T_h}
#' (default is \code{NULL} for both).
#' @param ch.data.pnts A logical argument for checking whether points \code{p1} and \code{p2} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num2PEtetra}}, \code{\link{Idom.num2PEtri}} and \code{\link{Idom.num2PEbasic.tri}}
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
#' n<-5 #try also n<-20
#' Xp<-runif.std.tetra(n)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))
#' r<-1.5
#'
#' Idom.num2PEstd.tetra(Xp[1,],Xp[2,],Xp,r)
#'
#' ind.gam2<-vector()
#' for (i in 1:(n-1))
#'  for (j in (i+1):n)
#'  {if (Idom.num2PEstd.tetra(Xp[i,],Xp[j,],Xp,r)==1)
#'   ind.gam2<-rbind(ind.gam2,c(i,j))}
#'
#' ind.gam2
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv;rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv
#' Idom.num2PEstd.tetra(Xp[1,],Xp[2,],Xp,r,rv1,rv2)
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv;
#' Idom.num2PEstd.tetra(Xp[1,],Xp[2,],Xp,r,rv1)
#'
#' #or try
#' rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv
#' Idom.num2PEstd.tetra(Xp[1,],Xp[2,],Xp,r,rv2=rv2)
#'
#' P1<-c(.1,.1,.1)
#' P2<-c(.4,.1,.2)
#' Idom.num2PEstd.tetra(P1,P2,Xp,r)
#'
#' Idom.num2PEstd.tetra(c(-1,-1,-1),Xp[2,],Xp,r,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE
#' #since not both points, p1 and p2, are data points in Xp
#' }
#'
#' @export
Idom.num2PEstd.tetra <- function(p1,p2,Xp,r,rv1=NULL,rv2=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1,3) || !is.point(p2,3))
  {stop('p1 and p2 must be a numeric 3D point')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('third argument must be numeric')}

  if (is.point(Xp,3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points are data points in Xp')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (isTRUE(all.equal(p1,p2)))
  {dom<-0; return(dom); stop}

  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  th<-rbind(A,B,C,D)  #standard regular tetrahedron

  if (is.null(rv1))
  {rv1<-rel.vert.tetraCC(p1,th)$rv #vertex region for point p1
  }
  if (is.null(rv2))
  {rv2<-rel.vert.tetraCC(p2,th)$rv  #vertex region for point p2
  }
  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IarcPEstd.tetra(p1,Xp[i,],r,rv1),IarcPEstd.tetra(p2,Xp[i,],r,rv2))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for three 3D points constituting a dominating set for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - standard regular tetrahedron case
#'
#' @description Returns \eqn{I(\{}\code{p1,p2,pt3}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp} in the standard regular
#' tetrahedron \eqn{T_h=T((0,0,0),(1,0,0),(1/2,\sqrt{3}/2,0),(1/2,\sqrt{3}/6,\sqrt{6}/3))}, that is,
#' returns 1 if \{\code{p1,p2,pt3}\} is a dominating set of PE-PCD, returns 0 otherwise.
#'
#' Point, \code{p1}, is in the region of vertex \code{rv1} (default is \code{NULL}), point, \code{p2}, is in the region of vertex \code{rv2}
#' (default is \code{NULL}); point, \code{pt3}), is in the region of vertex \code{rv3}) (default is \code{NULL}); vertices (and hence \code{rv1, rv2} and
#' \code{rv3}) are labeled as \code{1,2,3,4} in the order they are stacked row-wise in \eqn{T_h}.
#'
#' PE proximity region is constructed with
#' respect to the tetrahedron \eqn{T_h} with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center of mass \eqn{CM} (equivalent to circumcenter in this case).
#'
#'  \code{ch.data.pnts} is for checking whether points \code{p1}, \code{p2} and \code{pt3} are all data points in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the points \code{p1}, \code{p2} and \code{pt3} would constitute a dominating set
#' if they actually were all in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1,p2,pt3 Three 3D points to be tested for constituting a dominating set of the PE-PCD.
#' @param Xp A set of 3D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param rv1,rv2,rv3 The indices of the vertices whose regions contains \code{p1}, \code{p2} and \code{pt3},
#' respectively.
#' They take the vertex labels as \code{1,2,3,4} as in the row order of the vertices in \eqn{T_h}
#' (default is \code{NULL} for all).
#' @param ch.data.pnts A logical argument for checking whether points \code{p1} and \code{p2} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(\{}\code{p1,p2,pt3}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2,pt3}\} is a dominating set of PE-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num3PEtetra}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' tetra<-rbind(A,B,C,D)
#' n<-5 #try 20, 40, 100 (larger n may take a long time)
#' Xp<-runif.std.tetra(n)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))
#' r<-1.25
#'
#' Idom.num3PEstd.tetra(Xp[1,],Xp[2,],Xp[3,],Xp,r)
#'
#' ind.gam3<-vector()
#' for (i in 1:(n-2))
#'  for (j in (i+1):(n-1))
#'    for (k in (j+1):n)
#'  {if (Idom.num3PEstd.tetra(Xp[i,],Xp[j,],Xp[k,],Xp,r)==1)
#'   ind.gam3<-rbind(ind.gam3,c(i,j,k))}
#'
#' ind.gam3
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv; rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv;
#' rv3<-rel.vert.tetraCC(Xp[3,],tetra)$rv
#' Idom.num3PEstd.tetra(Xp[1,],Xp[2,],Xp[3,],Xp,r,rv1,rv2,rv3)
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv;
#' Idom.num3PEstd.tetra(Xp[1,],Xp[2,],Xp[3,],Xp,r,rv1)
#'
#' #or try
#' rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv
#' Idom.num3PEstd.tetra(Xp[1,],Xp[2,],Xp[3,],Xp,r,rv2=rv2)
#'
#' P1<-c(.1,.1,.1)
#' P2<-c(.3,.3,.3)
#' P3<-c(.4,.1,.2)
#' Idom.num3PEstd.tetra(P1,P2,P3,Xp,r)
#'
#' Idom.num3PEstd.tetra(Xp[1,],c(1,1,1),Xp[3,],Xp,r,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE since not all points are data points in Xp
#' }
#'
#' @export
Idom.num3PEstd.tetra <- function(p1,p2,pt3,Xp,r,rv1=NULL,rv2=NULL,rv3=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1,3) || !is.point(p2,3) || !is.point(pt3,3))
  {stop('p1, p2, and pt3 must be a numeric 3D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp,3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp) || !is.in.data(pt3,Xp))
    {stop('not all points are data points in Xp')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}


  if (isTRUE(all.equal(p1,p2)) || isTRUE(all.equal(p1,pt3)) || isTRUE(all.equal(p2,pt3)))
  {dom<-0; return(dom); stop}

  A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  th<-rbind(A,B,C,D)  #standard regular tetrahedron

  if (is.null(rv1))
  {rv1<-rel.vert.tetraCC(p1,th)$rv #vertex region for point p1
  }
  if (is.null(rv2))
  {rv2<-rel.vert.tetraCC(p2,th)$rv #vertex region for point p2
  }
  if (is.null(rv3))
  {rv3<-rel.vert.tetraCC(pt3,th)$rv  #vertex region for point pt3
  }

  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IarcPEstd.tetra(p1,Xp[i,],r,rv1),IarcPEstd.tetra(p2,Xp[i,],r,rv2),IarcPEstd.tetra(pt3,Xp[i,],r,rv3))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for a 3D point being a dominating point for Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' - one tetrahedron case
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 2D data set \code{Xp} in the tetrahedron \code{th}, that is,
#' returns 1 if \code{p} is a dominating point of PE-PCD, returns 0 otherwise.
#'
#' Point, \code{p}, is in the vertex region of vertex \code{rv} (default is \code{NULL}); vertices are labeled as \code{1,2,3,4}
#' in the order they are stacked row-wise in \code{th}.
#'
#' PE proximity region is constructed with respect to the tetrahedron \code{th} with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center of mass (\code{M="CM"}) or circumcenter (\code{M="CC"}) only.
#' and vertex regions are based on center of mass \eqn{CM} (equivalent to circumcenter in this case).
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p A 3D point that is to be tested for being a dominating point or not of the PE-PCD.
#' @param Xp A set of 3D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param rv Index of the vertex whose region contains point \code{p}, \code{rv} takes the vertex labels as \code{1,2,3,4} as
#' in the row order of the vertices in standard tetrahedron, default is \code{NULL}.
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num1PEstd.tetra}}, \code{\link{Idom.num1PEtri}} and \code{\link{Idom.num1PEbasic.tri}}
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
#' n<-5 #try also n<-20
#'
#' Xp<-runif.tetra(n,tetra)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))
#'
#' M<-"CM"; cent<-apply(tetra,2,mean)  #center of mass
#' #try also M<-"CC"; cent<-circumcenter.tetra(tetra)  #circumcenter
#'
#' r<-2
#'
#' P<-c(.4,.1,.2)
#' Idom.num1PEtetra(Xp[1,],Xp,tetra,r,M)
#' Idom.num1PEtetra(P,Xp,tetra,r,M)
#'
#' #or try
#' RV<-rel.vert.tetraCC(Xp[1,],tetra)$rv
#' Idom.num1PEtetra(Xp[1,],Xp,tetra,r,M,rv=RV)
#'
#' Idom.num1PEtetra(c(-1,-1,-1),Xp,tetra,r,M)
#' Idom.num1PEtetra(c(-1,-1,-1),c(-1,-1,-1),tetra,r,M)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1PEtetra(Xp[i,],Xp,tetra,r,M))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#' g1.pts<-Xp[ind.gam1,]
#'
#' Xlim<-range(tetra[,1],Xp[,1],cent[1])
#' Ylim<-range(tetra[,2],Xp[,2],cent[2])
#' Zlim<-range(tetra[,3],Xp[,3],cent[3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi =0,theta=40, bty = "g",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
#'          pch = 20, cex = 1, ticktype = "detailed")
#' #add the vertices of the tetrahedron
#' plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#' if (length(g1.pts)!=0)
#' {plot3D::points3D(g1.pts[,1],g1.pts[,2],g1.pts[,3], pch=4,col="red", add=TRUE)}
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#'
#' D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
#' L<-rbind(D1,D2,D3,D4,D5,D6); R<-rbind(cent,cent,cent,cent,cent,cent)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
#'
#' P<-c(.4,.1,.2)
#' Idom.num1PEtetra(P,Xp,tetra,r,M)
#'
#' Idom.num1PEtetra(c(-1,-1,-1),Xp,tetra,r,M,ch.data.pnt = FALSE)
#' #gives an error message if ch.data.pnt = TRUE since p is not a data point
#' }
#'
#' @export
Idom.num1PEtetra <- function(p,Xp,th,r,M="CM",rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p,3))
  {stop('p must be a numeric 3D point')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp,3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('point, p, is not a data point in Xp')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  if (isTRUE(all.equal(matrix(p,ncol=3),matrix(Xp,ncol=3))))
  {dom<-1;return(dom);stop}

  if (in.tetrahedron(p,th)$in.tetra==FALSE)
  {dom<-0; return(dom); stop}

  if (is.null(rv))
  {rv<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p,th)$rv,rel.vert.tetraCM(p,th)$rv); #vertex region for p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3,4))!=1)
  {stop('vertex index, rv, must be 1, 2, 3 or 4')}}

  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (IarcPEtetra(p,Xp[i,],th,r,M,rv=rv)==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for two 3D points constituting a dominating set for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - one tetrahedron case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp} in the tetrahedron \code{th}, that is,
#' returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD, returns 0 otherwise.
#'
#' Point, \code{p1}, is in the region of vertex \code{rv1} (default is \code{NULL}) and point, \code{p2}, is in the region of vertex \code{rv2}
#' (default is \code{NULL}); vertices (and hence \code{rv1} and \code{rv2}) are labeled as \code{1,2,3,4} in the order they are stacked
#' row-wise in \code{th}.
#'
#' PE proximity region is constructed with respect to the tetrahedron \code{th} with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center of mass (\code{M="CM"}) or circumcenter (\code{M="CC"}) only.
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1} and \code{p2} are both data points in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the points \code{p1} and \code{p2} would constitute a dominating set
#' if they actually were both in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1,p2 Two 3D points to be tested for constituting a dominating set of the PE-PCD.
#' @param Xp A set of 3D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param rv1,rv2 The indices of the vertices whose regions contains \code{p1} and \code{p2}, respectively.
#' They take the vertex labels as \code{1,2,3,4} as in the row order of the vertices in \code{th}
#' (default is \code{NULL} for both).
#' @param ch.data.pnts A logical argument for checking whether both points \code{p1} and \code{p2} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num2PEstd.tetra}}, \code{\link{Idom.num2PEtri}} and \code{\link{Idom.num2PEbasic.tri}}
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
#' n<-5
#'
#' set.seed(1)
#' Xp<-runif.tetra(n,tetra)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))
#'
#' M<-"CM";  #try also M<-"CC";
#' r<-1.5
#'
#' Idom.num2PEtetra(Xp[1,],Xp[2,],Xp,tetra,r,M)
#' Idom.num2PEtetra(c(-1,-1,-1),Xp[2,],Xp,tetra,r,M)
#'
#' ind.gam2<-ind.gamn2<-vector()
#' for (i in 1:(n-1))
#'  for (j in (i+1):n)
#'  {if (Idom.num2PEtetra(Xp[i,],Xp[j,],Xp,tetra,r,M)==1)
#'  {ind.gam2<-rbind(ind.gam2,c(i,j))
#'  }
#'  }
#' ind.gam2
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv;rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv
#' Idom.num2PEtetra(Xp[1,],Xp[2,],Xp,tetra,r,M,rv1,rv2)
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv;
#' Idom.num2PEtetra(Xp[1,],Xp[2,],Xp,tetra,r,M,rv1)
#'
#' #or try
#' rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv
#' Idom.num2PEtetra(Xp[1,],Xp[2,],Xp,tetra,r,M,rv2=rv2)
#'
#' P1<-c(.1,.1,.1)
#' P2<-c(.4,.1,.2)
#' Idom.num2PEtetra(P1,P2,Xp,tetra,r,M)
#'
#' Idom.num2PEtetra(c(-1,-1,-1),Xp[2,],Xp,tetra,r,M,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE
#' #since not both points, p1 and p2, are data points in Xp
#' }
#'
#' @export
Idom.num2PEtetra <- function(p1,p2,Xp,th,r,M="CM",rv1=NULL,rv2=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1,3) || !is.point(p2,3))
  {stop('p1 and p2 must be a numeric 3D point')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('third argument must be numeric')}

  if (is.point(Xp,3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points are data points in Xp')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  if (isTRUE(all.equal(p1,p2)))
  {dom<-0; return(dom); stop}

  if (is.null(rv1))
  {rv1<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p1,th)$rv,rel.vert.tetraCM(p1,th)$rv);} #vertex region for p1

  if (is.null(rv2))
  {rv2<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p2,th)$rv,rel.vert.tetraCM(p2,th)$rv);} #vertex region for p1

  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IarcPEtetra(p1,Xp[i,],th,r,M,rv1),IarcPEtetra(p2,Xp[i,],th,r,M,rv2))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for three 3D points constituting a dominating set for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - one tetrahedron case
#'
#' @description Returns \eqn{I(\{}\code{p1,p2,pt3}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp} in the tetrahedron \code{th}, that is,
#' returns 1 if \{\code{p1,p2,pt3}\} is a dominating set of PE-PCD, returns 0 otherwise.
#'
#' Point, \code{p1}, is in the region of vertex \code{rv1} (default is \code{NULL}), point, \code{p2}, is in the region of vertex \code{rv2}
#' (default is \code{NULL}); point, \code{pt3}), is in the region of vertex \code{rv3}) (default is \code{NULL}); vertices (and hence \code{rv1, rv2} and
#' \code{rv3}) are labeled as \code{1,2,3,4} in the order they are stacked row-wise in \code{th}.
#'
#' PE proximity region is constructed with
#' respect to the tetrahedron \code{th} with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center of mass \eqn{CM} (equivalent to circumcenter in this case).
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1}, \code{p2} and \code{pt3} are all data points in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the points \code{p1}, \code{p2} and \code{pt3} would constitute a dominating set
#' if they actually were all in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1,p2,pt3 Three 3D points to be tested for constituting a dominating set of the PE-PCD.
#' @param Xp A set of 3D points which constitutes the vertices of the PE-PCD.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#' @param rv1,rv2,rv3 The indices of the vertices whose regions contains \code{p1}, \code{p2} and \code{pt3},
#' respectively.
#' They take the vertex labels as \code{1,2,3,4} as in the row order of the vertices in \code{th}
#' ( default is \code{NULL} for all).
#' @param ch.data.pnts A logical argument for checking whether points \code{p1} and \code{p2} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(\{}\code{p1,p2,pt3}\} is a dominating set of the PE-PCD\eqn{)} where the vertices of the PE-PCD are the 3D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2,pt3}\} is a dominating set of PE-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num3PEstd.tetra}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' tetra<-rbind(A,B,C,D)
#' n<-5 #try 20, 40, 100 (larger n may take a long time)
#' Xp<-runif.tetra(n,tetra)$g
#'
#' M<-"CM";  #try also M<-"CC";
#' r<-1.25
#'
#' Idom.num3PEtetra(Xp[1,],Xp[2,],Xp[3,],Xp,tetra,r,M)
#'
#' ind.gam3<-vector()
#' for (i in 1:(n-2))
#'  for (j in (i+1):(n-1))
#'    for (k in (j+1):n)
#'    {if (Idom.num3PEtetra(Xp[i,],Xp[j,],Xp[k,],Xp,tetra,r,M)==1)
#'     ind.gam3<-rbind(ind.gam3,c(i,j,k))}
#'
#' ind.gam3
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv; rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv;
#' rv3<-rel.vert.tetraCC(Xp[3,],tetra)$rv
#' Idom.num3PEtetra(Xp[1,],Xp[2,],Xp[3,],Xp,tetra,r,M,rv1,rv2,rv3)
#'
#' #or try
#' rv1<-rel.vert.tetraCC(Xp[1,],tetra)$rv;
#' Idom.num3PEtetra(Xp[1,],Xp[2,],Xp[3,],Xp,tetra,r,M,rv1)
#'
#' #or try
#' rv2<-rel.vert.tetraCC(Xp[2,],tetra)$rv
#' Idom.num3PEtetra(Xp[1,],Xp[2,],Xp[3,],Xp,tetra,r,M,rv2=rv2)
#'
#' P1<-c(.1,.1,.1)
#' P2<-c(.3,.3,.3)
#' P3<-c(.4,.1,.2)
#' Idom.num3PEtetra(P1,P2,P3,Xp,tetra,r,M)
#'
#' Idom.num3PEtetra(Xp[1,],c(1,1,1),Xp[3,],Xp,tetra,r,M,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE since not all points are data points in Xp
#' }
#'
#' @export
Idom.num3PEtetra <- function(p1,p2,pt3,Xp,th,r,M="CM",rv1=NULL,rv2=NULL,rv3=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1,3) || !is.point(p2,3) || !is.point(pt3,3))
  {stop('p1, p2, and pt3 must be a numeric 3D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp,3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp) || !is.in.data(pt3,Xp))
    {stop('not all points are data points in Xp')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  if (isTRUE(all.equal(p1,p2)) || isTRUE(all.equal(p1,pt3)) || isTRUE(all.equal(p2,pt3)))
  {dom<-0; return(dom); stop}

  if (is.null(rv1))
  {rv1<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p1,th)$rv,rel.vert.tetraCM(p1,th)$rv);} #vertex region for p1

  if (is.null(rv2))
  {rv2<-ifelse(identical(M,"CC"),rel.vert.tetraCC(p2,th)$rv,rel.vert.tetraCM(p2,th)$rv);} #vertex region for p2

  if (is.null(rv3))
  {rv3<-ifelse(identical(M,"CC"),rel.vert.tetraCC(pt3,th)$rv,rel.vert.tetraCM(pt3,th)$rv);} #vertex region for pt3

  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IarcPEtetra(p1,Xp[i,],th,r,M,rv1),IarcPEtetra(p2,Xp[i,],th,r,M,rv2),
           IarcPEtetra(pt3,Xp[i,],th,r,M,rv3))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The domination number of Proportional Edge Proximity Catch Digraph (PE-PCD) - one tetrahedron case
#'
#' @description Returns the domination number of PE-PCD whose vertices are the data points in \code{Xp}.
#'
#' PE proximity region is defined with respect to the tetrahedron \code{th} with expansion parameter \eqn{r \ge 1} and
#' vertex regions are based on the center \code{M} which is circumcenter (\code{"CC"}) or center of mass (\code{"CM"}) of \code{th}
#' with default=\code{"CM"}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 3D points which constitute the vertices of the digraph.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param r A positive real number which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#'
#' @return A \code{list} with two elements
#' \item{dom.num}{Domination number of PE-PCD with vertex set = \code{Xp} and expansion parameter \eqn{r \ge 1} and center \code{M}}
#' \item{mds}{A minimum dominating set of PE-PCD with vertex set = \code{Xp} and expansion parameter \eqn{r \ge 1} and center \code{M}}
#' \item{ind.mds}{Indices of the minimum dominating set \code{mds}}
#'
#' @seealso \code{\link{PEdom.num.tri}}
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
#' n<-10  #try also n<-20
#'
#' Xp<-runif.tetra(n,tetra)$g
#'
#' M<-"CM"  #try also M<-"CC"
#' r<-1.25
#'
#' PEdom.num.tetra(Xp,tetra,r,M)
#'
#' P1<-c(.5,.5,.5)
#' PEdom.num.tetra(P1,tetra,r,M)
#' }
#'
#' @export PEdom.num.tetra
PEdom.num.tetra <- function(Xp,th,r,M="CM")
{
  if (!is.numeric(as.matrix(Xp)) )
  {stop('Xp must be numeric')}

  if (is.point(Xp,dim=3))
  { Xp<-matrix(Xp,ncol=3)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=3 )
  {stop('Xp must be of dimension nx3')}
  }

  th<-as.matrix(th)
  if (!is.numeric(th) || nrow(th)!=4 || ncol(th)!=3)
  {stop('th must be numeric and of dimension 4x3')}

  vec1<-rep(1,4);
  D0<-det(matrix(cbind(th,vec1),ncol=4))
  if (round(D0,14)==0)
  {stop('the tetrahedron is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  n<-nrow(Xp)  #number of Xp points

  ind.tetra<-mds<-mds.ind<-c()
  for (i in 1:n)
  {
    if (in.tetrahedron(Xp[i,],th,boundary = TRUE)$i)
      ind.tetra<-c(ind.tetra,i)  #indices of data points in the tetrahedron, th
  }

  Xth<-matrix(Xp[ind.tetra,],ncol=3) #data points in the tetrahedron, th

  nth<-nrow(Xth)  #number of points inside the tetrahedron
  if (nth==0)
  {gam<-0;
  res<-list(dom.num=gam, #domination number
            mds=NULL,  #a minimum dominating set
            ind.mds=NULL #indices of the mds
  )
  return(res); stop}

  Cl2f0<-cl2faces.vert.reg.tetra(Xth,th,M)
  Cl2f<-Cl2f0$ext #for general r, points closest to opposite edges in the vertex regions
  Cl2f.ind<-Cl2f0$ind # indices of these extrema wrt Xth
  Ext.ind =ind.tetra[Cl2f.ind] #indices of these extrema wrt to the original data

  #Gamma=1 piece
  cnt<-0; j<-1;
  while (j<=4 & cnt==0)
  {
    if (sum(!is.na(Cl2f[j,]))==0 )
    {j<-j+1
    } else
    {
      if (Idom.num1PEtetra(Cl2f[j,],Xth,th,r,M,rv=j)==1)
      {gam<-1; cnt<-1; mds<-rbind(mds,Cl2f[j,]); mds.ind=c(mds.ind,Ext.ind[j])
      } else
      {j<-j+1}
    }
  }

  #Gamma=2 piece
  if (cnt==0)
  { k<-1; cnt2<-0;
  while (k<=3 & cnt2==0)
  {l<-k+1;
  while (l<=4 & cnt2==0)
  {
    if (sum(!is.na(Cl2f[k,]))==0 | sum(!is.na(Cl2f[l,]))==0 )
    {l<-l+1
    } else
    {
      if (Idom.num2PEtetra(Cl2f[k,],Cl2f[l,],Xth,th,r,M,rv1=k,rv2=l)==1)
      {gam<-2;cnt2<-1; mds<-rbind(mds,Cl2f[c(k,l),]); mds.ind=c(mds.ind,Ext.ind[c(k,l)])
      } else {l<-l+1};
    }
  }
  k<-k+1;
  }
  }

  #Gamma=3 piece
  if (cnt==0 && cnt2==0)
  { k3<-1; cnt3<-0;
  while (k3<=2 & cnt3==0)
  {l3<-k3+1;
  while (l3<=3 & cnt3==0)
  {m3<-l3+1;
  while (m3<=4 & cnt3==0)
  {
    if (sum(!is.na(Cl2f[k3,]))==0 | sum(!is.na(Cl2f[l3,]))==0 | sum(!is.na(Cl2f[m3,]))==0 )
    {m3<-m3+1
    } else
    {
      if (Idom.num3PEtetra(Cl2f[k3,],Cl2f[l3,],Cl2f[m3,],Xth,th,r,M,rv1=k3,rv2=l3,rv3=m3)==1)
      {gam<-3; cnt3<-1; mds<-rbind(mds,Cl2f[c(k3,l3,m3),]); mds.ind=c(mds.ind,Ext.ind[c(k3,l3,m3)])
      } else {m3<-m3+1};
    }
  }
  l3<-l3+1;
  }
  k3<-k3+1
  }
  }

  if (cnt==0 && cnt2==0 && cnt3==0)
  {gam <-4; mds<-rbind(mds,Cl2f); mds.ind=c(mds.ind,Ext.ind)}

  row.names(mds)<-c()
  res<-list(dom.num=gam, #domination number
       mds=mds, #a minimum dominating set
       ind.mds =mds.ind #indices of a minimum dominating set (wrt to original data)
  )
res
} #end of the function
#'
