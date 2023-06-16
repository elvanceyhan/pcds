#CentSim2D.R;
#Functions for NCS in R^2
#################################################################

#Contains the functions used in PCD calculations, such as generation of data in a given triangle
#estimation of gamma, arc density etc.

# funsMuVarCS2D
#'
#' @title Returns the mean and (asymptotic) variance of arc density of  Central Similarity Proximity Catch Digraph (CS-PCD)
#' for 2D uniform data in one triangle
#'
#' @description
#' Two functions: \code{muCS2D} and \code{asyvarCS2D}.
#'
#' \code{muCS2D} returns the mean of the (arc) density of CS-PCD
#' and \code{asyvarCS2D} returns the asymptotic variance of the arc density of CS-PCD
#' with expansion parameter \eqn{t>0} for 2D uniform data in a triangle.
#'
#' CS proximity regions are defined with respect to the triangle and
#' vertex regions are based on center of mass, \eqn{CM} of the triangle.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#'
#' @return \code{muCS2D} returns the mean and \code{asyvarCS2D} returns the (asymptotic) variance of the
#' arc density of CS-PCD for uniform data in any triangle
#'
#' @name funsMuVarCS2D
NULL
#'
#' @seealso \code{\link{muPE2D}} and \code{\link{asyvarPE2D}}
#'
#' @rdname funsMuVarCS2D
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for muCS2D
#' muCS2D(.5)
#'
#' tseq<-seq(0.01,5,by=.1)
#' ltseq<-length(tseq)
#'
#' mu<-vector()
#' for (i in 1:ltseq)
#' {
#'   mu<-c(mu,muCS2D(tseq[i]))
#' }
#'
#' plot(tseq, mu,type="l",xlab="t",ylab=expression(mu(t)),lty=1,xlim=range(tseq))
#' }
#'
#' @export
muCS2D <- function(t)
{
  if (!is.point(t,1) || t<=0)
  {stop('the argument must be a scalar greater than 0')}

  mean<-0;
  if (t < 1)
  {
    mean<-t^2/6;
  } else {
    mean<-(t*(4*t-1))/(2*(1+2*t)*(2+t));
  }
  mean
} #end of the function
#'
#' @rdname funsMuVarCS2D
#'
#' @examples
#' \dontrun{
#' #Examples for asyvarCS2D
#' asyvarCS2D(.5)
#'
#' tseq<-seq(0.01,10,by=.1)
#' ltseq<-length(tseq)
#'
#' asyvar<-vector()
#' for (i in 1:ltseq)
#' {
#'   asyvar<-c(asyvar,asyvarCS2D(tseq[i]))
#' }
#'
#' par(mar=c(5,5,4,2))
#' plot(tseq, asyvar,type="l",xlab="t",ylab=expression(paste(sigma^2,"(t)")),lty=1,xlim=range(tseq))
#' }
#'
#' @export
asyvarCS2D <- function(t)
{
  if (!is.point(t,1) || t<=0)
  {stop('the argument must be a scalar greater than 0')}

  asyvar<-0;
  if (t < 1)
  {
    asyvar<-(t^4*(6*t^5-3*t^4-25*t^3+t^2+49*t+14))/(45*(t+1)*(2*t+1)*(t+2)) ;
  } else {
    asyvar<-(168*t^7+886*t^6+1122*t^5+45*t^4-470*t^3-114*t^2+48*t+16)/(5*(2*t+1)^4*(t+2)^4);
  }
  asyvar
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Central Similarity Proximity Catch
#' Digraphs (CS-PCDs) - first one-sixth of the standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1))} for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t=1)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t=1)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{t=1}.
#'
#' CS proximity region is defined with respect to the standard equilateral triangle
#' \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} and edge regions are based on the center of mass \eqn{CM=(1/2,\sqrt{3}/6)}.
#' Here \code{p1} must lie in the first one-sixth of \eqn{T_e}, which is the triangle with vertices \eqn{T(A,D_3,CM)=T((0,0),(1/2,0),CM)}.
#' If \code{p1} and \code{p2} are distinct and \code{p1} is outside of \eqn{T(A,D_3,CM)} or \code{p2} is outside \eqn{T_e}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1))} for \code{p1} in the first one-sixth of \eqn{T_e},
#' \eqn{T(A,D_3,CM)}, that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t=1)}, returns 0 otherwise
#'
#' @seealso \code{\link{IarcCSstd.tri}}
#'
#' @author Elvan Ceyhan
#'
IarcCS.Te.onesixth <- function(p1,p2)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('both arguments must be numeric 2D points')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  CM<-(A+B+C)/3; D3<-(A+B)/2;
  tri<-rbind(A,D3,CM)
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,rbind(A,B,C), boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  arc<-0
  Dx<--1.732050807*y[1]+x[1]; Ex<-1.732050807*y[1]+x[1];
  if (y[2] < 1.732050808*x[2]-1.732050808*Dx && y[2] < -1.732050808*x[2]+1.732050808*Ex)
  {
    arc<-1
  }
  arc
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' - first one-sixth of the standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of the 2D data set \code{Xp} of CS-PCD\eqn{)} in the standard equilateral
#' triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}, that is, returns 1 if \code{p} is a dominating point of CS-PCD,
#' returns 0 otherwise.
#'
#' Point, \code{p}, must lie in the first one-sixth of \eqn{T_e}, which is the triangle with vertices
#' \eqn{T(A,D_3,CM)=T((0,0),(1/2,0),CM)}.
#'
#' CS proximity region is constructed with respect to \eqn{T_e} with expansion parameter \eqn{t=1}.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks whether the point \code{p} would be a dominating point if it actually were in the data
#' set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point or not of the CS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the CS-PCD.
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num1CSstd.tri}} and \code{\link{Idom.num1CSt1std.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
Idom.num1CS.Te.onesixth <- function(p,Xp,ch.data.pnt=FALSE)
{
  if (!is.point(p))
  {stop('p must be a numeric point of dimension 2')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('point, p, is not a data point in Xp')}
  }

  n<-nrow(Xp)
  ct<-1; i<-1;
  while (i <= n && ct==1)
  {if (IarcCS.Te.onesixth(p,Xp[i,])==0)
  {ct<-0};
    i<-i+1;
  }
  ct
} #end of the function
#'

#################################################################

#' @title The indicator for two points constituting a dominating set for Central Similarity Proximity Catch Digraphs
#' (CS-PCDs) - first one-sixth of the standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is a dominating set of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp}),
#' that is, returns 1 if \code{p} is a dominating point of CS-PCD, returns 0 otherwise.
#'
#' CS proximity region is
#' constructed with respect to the standard equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} and
#' with expansion parameter \eqn{t=1}. Point, \code{p1}, must lie in the first one-sixth of \eqn{T_e}, which is the triangle with
#' vertices \eqn{T(A,D_3,CM)=T((0,0),(1/2,0),CM)}.
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1} and \code{p2} are data points in \code{Xp} or not
#' (default is \code{FALSE}), so by default this function checks whether the points \code{p1} and \code{p2} would be a
#' dominating set if they actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis;textual}{pcds}).
#'
#' @param p1,p2 Two 2D points to be tested for constituting a dominating set of the CS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the CS-PCD.
#' @param ch.data.pnts A logical argument for checking whether points \code{p1} and \code{p2} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp}),
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of CS-PCD, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num2CSstd.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
Idom.num2CS.Te.onesixth <- function(p1,p2,Xp,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points are data points in Xp')}
  }

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  CM<-(A+B+C)/3; D3<-(A+B)/2;
  tri<-rbind(A,D3,CM); Te<-rbind(A,B,C)

  n<-nrow(Xp)
  ct<-1; i<-1;
  while (i<=n && ct==1)
  {
    if (max(IarcCS.Te.onesixth(p1,Xp[i,]),IarcCS.Te.onesixth(p2,Xp[i,]))==0)
    {ct<-0};
    i<-i+1;
  }
  ct
} #end of the function
#'

#################################################################

#' @title The vertices of the Central Similarity (CS) Proximity Region in a general triangle
#'
#' @description Returns the vertices of the CS proximity region (which is itself a triangle) for a point in the
#' triangle \code{tri}\eqn{=T(A,B,C)=}\code{(rv=1,rv=2,rv=3)}.
#'
#' CS proximity region is defined with respect to the triangle \code{tri}
#' with expansion parameter \eqn{t>0} and edge regions based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' Edge regions are labeled as \code{1,2,3} rowwise for the corresponding vertices
#' of the triangle \code{tri}. \code{re} is the index of the edge region \code{p} resides, with default=\code{NULL}.
#' If \code{p} is outside of \code{tri}, it returns \code{NULL} for the proximity region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p A 2D point whose CS proximity region is to be computed.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' @param re Index of the \code{M}-edge region containing the point \code{p},
#' either \code{1,2,3} or \code{NULL} (default is \code{NULL}).
#'
#' @return Vertices of the triangular region which constitutes the CS proximity region with expansion parameter
#' \eqn{t>0} and center \code{M} for a point \code{p}
#'
#' @seealso \code{\link{NPEtri}}, \code{\link{NAStri}}, and  \code{\link{IarcCStri}}
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
#' tau<-1.5
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' n<-3
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' NCStri(Xp[1,],Tr,tau,M)
#'
#' P1<-as.numeric(runif.tri(1,Tr)$g)  #try also P1<-c(.4,.2)
#' NCStri(P1,Tr,tau,M)
#'
#' #or try
#' re<-rel.edges.tri(P1,Tr,M)$re
#' NCStri(P1,Tr,tau,M,re)
#' }
#'
#' @export NCStri
NCStri <- function(p,tri,t,M=c(1,1,1),re=NULL)
{
  if (!is.point(p) )
  {stop('p must be a numeric 2D point')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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

  if (!in.triangle(p,tri,boundary=TRUE)$in.tri)
  {reg<-NULL; return(reg); stop}

  if (is.null(re))
  {re<-rel.edge.tri(p,tri,M)$re #edge region for p
  } else
  {  if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')}}

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,];

  #If p is on the boundary of the triangle
  if ( min(dist.point2line(p,A,B)$dist,dist.point2line(p,A,C)$dist,dist.point2line(p,B,C)$dist)==0)
  {reg<-rbind(p,p,p); row.names(reg)<-c(); return(reg); stop}

  if (re==1)
  {d1<-t*dist.point2line(p,B,C)$d; d2<-dist.point2line(M,B,C)$d;
  sr<-d1/d2;
  tri.shr<-sr*(tri-t(replicate(3,M)))
  reg<-tri.shr+t(replicate(3,p))
  A1<-reg[1,]; B1<-reg[2,]; C1<-reg[3,];
  if (t>1)
  {G1<-intersect2lines(A1,B1,B,C);
  H1<-intersect2lines(A1,C1,B,C);
  if (in.triangle(A1,tri,boundary = TRUE)$i)
  {reg<-rbind(A1,G1,H1)
  } else if (in.triangle(C,reg,boundary = FALSE)$i && !in.triangle(B,reg,boundary = TRUE)$i)
  {
    I1<-intersect2lines(A1,B1,A,C);
    reg<-rbind(G1,C,I1)
  } else if (!in.triangle(C,reg,boundary = FALSE)$i && in.triangle(B,reg,boundary = TRUE)$i)
  {J1<-intersect2lines(A1,C1,A,B)
  reg<-rbind(H1,B,J1)
  } else
  {reg<-tri}
  }
  } else if (re==2)
  {
    d1<-t*dist.point2line(p,A,C)$d; d2<-dist.point2line(M,A,C)$d;
    sr<-d1/d2;
    tri.shr<-sr*(tri-t(replicate(3,M)))
    reg<-tri.shr+t(replicate(3,p))
    A1<-reg[1,]; B1<-reg[2,]; C1<-reg[3,];
    if (t>1)
    {G1<-intersect2lines(A1,B1,A,C);
    H1<-intersect2lines(B1,C1,A,C);
    if (in.triangle(B1,tri,boundary = TRUE)$i)
    {reg<-rbind(B1,G1,H1)
    } else if (in.triangle(C,reg,boundary = FALSE)$i && !in.triangle(A,reg,boundary = TRUE)$i)
    {
      I1<-intersect2lines(A1,B1,B,C);
      reg<-rbind(G1,C,I1)
    } else if (!in.triangle(C,reg,boundary = FALSE)$i && in.triangle(A,reg,boundary = TRUE)$i)
    {J1<-intersect2lines(B1,C1,A,B)
    reg<-rbind(H1,A,J1)
    } else
    {reg<-tri}
    }
  } else
  {
    d1<-t*dist.point2line(p,A,B)$d; d2<-dist.point2line(M,A,B)$d;
    sr<-d1/d2;
    tri.shr<-sr*(tri-t(replicate(3,M)))
    reg<-tri.shr+t(replicate(3,p))
    A1<-reg[1,]; B1<-reg[2,]; C1<-reg[3,];
    if (t>1)
    {G1<-intersect2lines(A1,C1,A,B);
    H1<-intersect2lines(B1,C1,A,B);
    if (in.triangle(C1,tri,boundary = TRUE)$i)
    {reg<-rbind(C1,G1,H1)
    } else if (in.triangle(B,reg,boundary = FALSE)$i && !in.triangle(A,reg,boundary = TRUE)$i)
    {
      I1<-intersect2lines(A1,C1,B,C);
      reg<-rbind(G1,B,I1)
    } else if (!in.triangle(B,reg,boundary = FALSE)$i && in.triangle(A,reg,boundary = TRUE)$i)
    {J1<-intersect2lines(B1,C1,A,C)
    reg<-rbind(H1,A,J1)
    } else
    {reg<-tri}
    }
  }
  row.names(reg)<-c()
  reg
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from one point to another for Central Similarity Proximity
#' Catch Digraphs (CS-PCDs)
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2}, that is,
#' returns 1 if \code{p2} is in \eqn{NCS(p1,t)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with the expansion parameter \eqn{t>0}.
#'
#' CS proximity region is constructed with respect to the triangle \code{tri} and
#' edge regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \code{tri}
#' or based on the circumcenter of \code{tri}.
#' \code{re} is the index of the edge region \code{p} resides, with default=\code{NULL}
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \code{tri}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}.
#' @param re Index of the \code{M}-edge region containing the point \code{p},
#' either \code{1,2,3} or \code{NULL} (default is \code{NULL}).
#'
#' @return I(\code{p2} is in \eqn{NCS(p1,t)}) for \code{p1}, that is, returns 1 if \code{p2} is in \eqn{NCS(p1,t)}, returns 0 otherwise
#'
#' @seealso \code{\link{IarcAStri}}, \code{\link{IarcPEtri}}, \code{\link{IarcCStri}}, and \code{\link{IarcCSstd.tri}}
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
#' tau<-1.5
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.2)
#'
#' n<-10
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' IarcCStri(Xp[1,],Xp[2,],Tr,tau,M)
#'
#' P1<-as.numeric(runif.tri(1,Tr)$g)
#' P2<-as.numeric(runif.tri(1,Tr)$g)
#' IarcCStri(P1,P2,Tr,tau,M)
#'
#' #or try
#' re<-rel.edges.tri(P1,Tr,M)$re
#' IarcCStri(P1,P2,Tr,tau,M,re)
#' }
#'
#' @export IarcCStri
IarcCStri <- function(p1,p2,tri,t,M,re=NULL)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,tri,boundary=TRUE)$in.tri || !in.triangle(p2,tri,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(re))
  {rel.ed<-rel.edge.tri(p1,tri,M)$re #rel.edges.tri(p1,tri,M)$re #edge region for p1
  } else if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')
  } else
  {rel.ed<-re}

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,];

  if (rel.ed==1)
  {d1<-t*dist.point2line(p1,B,C)$d; d2<-dist.point2line(M,B,C)$d;
  sr<-d1/d2;
  tri.shr<-sr*(tri-t(replicate(3,M)))
  reg<-tri.shr+t(replicate(3,p1))
  } else if (rel.ed==2)
  {
    d1<-t*dist.point2line(p1,A,C)$d; d2<-dist.point2line(M,A,C)$d;
    sr<-d1/d2;
    tri.shr<-sr*(tri-t(replicate(3,M)))
    reg<-tri.shr+t(replicate(3,p1))
  } else
  {
    d1<-t*dist.point2line(p1,A,B)$d; d2<-dist.point2line(M,A,B)$d;
    sr<-d1/d2;
    tri.shr<-sr*(tri-t(replicate(3,M)))
    reg<-tri.shr+t(replicate(3,p1))
  }

  arc<-ifelse(area.polygon(reg)==0,0,sum(in.triangle(p2,reg,boundary=TRUE)$in.tri))
  arc
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Central Similarity Proximity Catch
#' Digraphs (CS-PCDs) - standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{r \ge 1}.
#'
#' CS proximity region is defined with respect to the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Edge regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
#' barycentric coordinates in the interior of the standard basic triangle \eqn{T_b};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
#' \code{re} is the index of the edge region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_b}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions (i.e., translation, rotation, and reflection) and scaling,
#' preserving uniformity of the points in the original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region; must be \eqn{\ge 1}
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle or circumcenter of \eqn{T_b};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
#' @param re The index of the edge region in \eqn{T_b} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise
#'
#' @seealso \code{\link{IarcCStri}} and \code{\link{IarcCSstd.tri}}
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
#' M<-as.numeric(runif.basic.tri(1,c1,c2)$g)
#'
#' tau<-2
#'
#' P1<-as.numeric(runif.basic.tri(1,c1,c2)$g)
#' P2<-as.numeric(runif.basic.tri(1,c1,c2)$g)
#' IarcCSbasic.tri(P1,P2,tau,c1,c2,M)
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' IarcCSbasic.tri(P1,P2,tau,c1,c2,M)
#' IarcCSbasic.tri(P1,P1,tau,c1,c2,M)
#'
#' #or try
#' Re<-rel.edge.basic.tri(P1,c1,c2,M)$re
#' IarcCSbasic.tri(P1,P2,tau,c1,c2,M,Re)
#' IarcCSbasic.tri(P1,P1,tau,c1,c2,M,Re)
#' }
#'
#' @export IarcCSbasic.tri
IarcCSbasic.tri <- function(p1,p2,t,c1,c2,M=c(1,1,1),re=NULL)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

  if (!is.point(t,1) || t<0)
  {stop('t must be a scalar > 0')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1 < 0 || c1 > 1/2 || c2 <= 0 || (1-c1)^2+c2^2 > 1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates
          or the circumcenter "CC" ')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,Tb,boundary=TRUE)$in.tri || !in.triangle(p2,Tb,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(re))
  { re<-rel.edge.basic.tri(p1,c1,c2,M)$re
  } else
  {  if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')}}

  X1<-p1[1]; Y1<-p1[2];
  X2<-p2[1]; Y2<-p2[2];
  m1<-M[1]; m2<-M[2];

  if (re==1)
  {d1<-t*abs(c2*X1+(1-c1)*Y1-c2)/sqrt((1-c1)^2+c2^2); d2<-abs(c2*m1+(1-c1)*m2-c2)/sqrt((1-c1)^2+c2^2);
  sr<-d1/d2;
  tri.shr<-sr*(Tb-t(replicate(3,M)))
  reg<-tri.shr+t(replicate(3,p1))
  } else if (re==2)
  {
    d1<-t*abs(c2*X1-c1*Y1)/sqrt(c1^2+c2^2); d2<-abs(c2*m1-c1*m1)/sqrt(c1^2+c2^2);
    sr<-d1/d2;
    tri.shr<-sr*(Tb-t(replicate(3,M)))
    reg<-tri.shr+t(replicate(3,p1))
  } else
  {
    d1<-t*Y1; d2<-m2;
    sr<-d1/d2;
    tri.shr<-sr*(Tb-t(replicate(3,M)))
    reg<-tri.shr+t(replicate(3,p1))
  }

  arc<-ifelse(area.polygon(reg)==0,0,sum(in.triangle(p2,reg,boundary=TRUE)$in.tri))
  arc
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Central Similarity Proximity Catch
#' Digraphs (CS-PCDs) - standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{t >0}.
#'
#' CS proximity region is defined with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e}.
#' \code{rv} is the index of the vertex region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_e}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param re The index of the edge region in \eqn{T_e} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise
#'
#' @seealso \code{\link{IarcCStri}}, \code{\link{IarcCSbasic.tri}}, and \code{\link{IarcPEstd.tri}}
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
#' n<-3
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2) or M=(A+B+C)/3
#'
#' IarcCSstd.tri(Xp[1,],Xp[3,],t=2,M)
#' IarcCSstd.tri(c(0,1),Xp[3,],t=2,M)
#'
#' #or try
#' Re<-rel.edge.tri(Xp[1,],Te,M) $re
#' IarcCSstd.tri(Xp[1,],Xp[3,],t=2,M,Re)
#' }
#'
#' @export IarcCSstd.tri
IarcCSstd.tri <- function(p1,p2,t,M=c(1,1,1),re=NULL)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

  if (!is.point(t,1) || t<0)
  {stop('t must be a scalar >0')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}
  m1=M[1]; m2=M[2]

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,Te,boundary=TRUE)$in.tri || !in.triangle(p2,Te,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(re))
  {re<-rel.edge.tri(p1,Te,M)$re
  } else
  {  if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')}}

  x1<-p1[1]; y1<-p1[2];
  # x2<-p2[1]; y2<-p2[2];
  arc<-0;
  if (re==1)
  {
    D1x=(-m1*t*y1 + m2*x1)/m2; D1y= -t*y1 + y1
    E1x=(m2*x1 - t*y1*(m1 - 1))/m2; E1y=-t*y1 + y1
    F1x=(2*m2*x1 - 2*t*y1*(m1 - 1/2))/(2*m2); F1y=(t*sqrt(3)*y1 - 2*y1*m2*(t - 1))/(2*m2)
    reg=matrix(c(D1x,D1y,E1x,E1y,F1x,F1y),ncol=2,byrow = TRUE)
  } else {
    if (re==2)
    {
      D1x=x1 + t*((m1*sqrt(3) - y1*m1 + m2*x1 - sqrt(3) + y1)/(m1*sqrt(3) - sqrt(3) + m2) - x1);
      D1y= y1 + t*(sqrt(3)*((1 - x1)*m2 + y1*(m1 - 1))/((m1 - 1)*sqrt(3) + m2) - y1)
      E1x=x1 + t*((2*m1*sqrt(3) - sqrt(3)*x1 - 2*y1*m1 + 2*m2*x1 - sqrt(3) + y1)/(2*(m1*sqrt(3) - sqrt(3) + m2)) - x1);
      E1y=y1 + t*(sqrt(3)*(sqrt(3)*x1 + 2*y1*m1 - 2*m2*x1 - sqrt(3) + 2*m2 - y1)/((2*m1 - 2)*sqrt(3) + 2*m2) - y1)
      F1x=-(sqrt(3)*m1*t*y1 - sqrt(3)*m2*x1 + 3*m1*t*x1 - 3*m1*t - 3*m1*x1 + 3*x1)/(sqrt(3)*m2 + 3*m1 - 3);
      F1y=y1 + t*(sqrt(3)*((1 - x1)*m2 + y1*(m1 - 1))/((m1 - 1)*sqrt(3) + m2) - y1)
      reg=matrix(c(D1x,D1y,E1x,E1y,F1x,F1y),ncol=2,byrow = TRUE)
    } else {
      D1x=(-x1*m1*(t - 1)*sqrt(3) + t*y1*m1 - m2*x1)/(m1*sqrt(3) - m2);
      D1y= ((-m2*t*x1 + m1*y1)*sqrt(3) + y1*m2*(t - 1))/(m1*sqrt(3) - m2)
      E1x=(-2*x1*(m1*(t - 1) - t/2)*sqrt(3) + 2*t*y1*m1 - t*y1 - 2*m2*x1)/(2*m1*sqrt(3) - 2*m2);
      E1y=((-2*m2*t*x1 - y1*(t - 2*m1))*sqrt(3) + 2*y1*m2*(t - 1) + 3*t*x1)/(2*m1*sqrt(3) - 2*m2)
      F1x=((m1*t*y1 - m2*x1 - t*y1)*sqrt(3) - 3*x1*(m1*(t - 1) - t))/(-sqrt(3)*m2 + 3*m1);
      F1y=((-m2*t*x1 + m1*y1)*sqrt(3) + y1*m2*(t - 1))/(m1*sqrt(3) - m2)
      reg=matrix(c(D1x,D1y,E1x,E1y,F1x,F1y),ncol=2,byrow = TRUE)
    }}

  arc<-ifelse(area.polygon(reg)==0,0,sum(in.triangle(p2,reg,boundary=TRUE)$in.tri))
  arc
} #end of the function
#'

#################################################################

# funsCSt1EdgeRegs
#'
#' @title Each function is for the presence of an arc from a point in one of the edge regions
#' to another for Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' standard equilateral triangle case with \eqn{t=1}
#'
#' @description
#' Three indicator functions: \code{IarcCSt1.std.triRAB}, \code{IarcCSt1.std.triRBC} and \code{IarcCSt1.std.triRAC}.
#'
#' The function \code{IarcCSt1.std.triRAB} returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1)} for \code{p1} in \eqn{RAB}
#' (edge region for edge \eqn{AB}, i.e., edge 3) in the standard equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))};
#'
#' \code{IarcCSt1.std.triRBC} returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1)} for \code{p1} in \eqn{RBC} (edge region for edge \eqn{BC}, i.e., edge 1) in \eqn{T_e};
#' and
#'
#' \code{IarcCSt1.std.triRAC} returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1)} for \code{p1} in \eqn{RAC} (edge region for edge \eqn{AC}, i.e., edge 2) in \eqn{T_e}.
#'
#' That is, each function returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t=1)}, returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the
#' CS proximity region for point \eqn{x} with expansion parameter \eqn{t=1}.
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#'
#' @return Each function returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1))} for \code{p1}, that is,
#' returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t=1)}, returns 0 otherwise
#'
#' @name funsCSt1EdgeRegs
NULL
#'
#' @seealso \code{\link{IarcCSstd.triRAB}}, \code{\link{IarcCSstd.triRBC}} and \code{\link{IarcCSstd.triRAC}}
#'
#' @rdname funsCSt1EdgeRegs
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for IarcCSt1.std.triRAB
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' T3<-rbind(A,B,CM);
#'
#' set.seed(1)
#' Xp<-runif.std.tri(10)$gen.points
#'
#' IarcCSt1.std.triRAB(Xp[1,],Xp[2,])
#'
#' IarcCSt1.std.triRAB(c(.2,.5),Xp[2,])
#' }
#'
#' @export
IarcCSt1.std.triRAB <- function(p1,p2)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('both arguments must be numeric 2D points')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  CM<-(A+B+C)/3;
  tri<-rbind(A,B,CM); Te<-rbind(A,B,C)
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  arc<-0
  if (y[2] < 1.732050808*x[2] + 3.0*y[1] - 1.732050808*x[1]
      &&
      y[2] < -1.732050808*x[2] + 3.0*y[1] + 1.732050808*x[1])
  {
    arc<-1
  }
  arc
} #end of the function
#'
#' @rdname funsCSt1EdgeRegs
#'
#' @examples
#' \dontrun{
#' #Examples for IarcCSt1.std.triRBC
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' T1<-rbind(B,C,CM);
#'
#' set.seed(1)
#' Xp<-runif.std.tri(3)$gen.points
#'
#' IarcCSt1.std.triRBC(Xp[1,],Xp[2,])
#'
#' IarcCSt1.std.triRBC(c(.2,.5),Xp[2,])
#' }
#'
#' @export
IarcCSt1.std.triRBC <- function(p1,p2)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('both arguments must be numeric 2D points')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  CM<-(A+B+C)/3;
  tri<-rbind(B,C,CM); Te<-rbind(A,B,C)
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  arc<-0;
  Dy <--.8660254040+1.500000000*y[1]+.8660254042*x[1];
  if (y[2] > Dy && y[2]< 1.732050808*x[2]+1.732050808-3.464101616*x[1])
  {arc<-1};
  arc
} #end of the function
#'
#' @rdname funsCSt1EdgeRegs
#'
#' @examples
#' \dontrun{
#' #Examples for IarcCSt1.std.triRAC
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' T2<-rbind(A,C,CM);
#'
#' set.seed(1)
#' Xp<-runif.std.tri(3)$gen.points
#'
#' IarcCSt1.std.triRAC(Xp[1,],Xp[2,])
#' IarcCSt1.std.triRAC(c(1,2),Xp[2,])
#' }
#'
#' @export
IarcCSt1.std.triRAC <- function(p1,p2)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('both arguments must be numeric 2D points')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  CM<-(A+B+C)/3;
  tri<-rbind(A,C,CM); Te<-rbind(A,B,C)
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  arc<-0;
  Dy <-1.500000000*y[1]-.8660254042*x[1];
  if (y[2] > Dy && y[2]< -1.732050808*x[2]+3.464101616*x[1])
  {arc<-1};
  arc
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Central Similarity Proximity Catch
#' Digraphs (CS-PCDs) - standard equilateral triangle case with \eqn{t=1}
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1))} for points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t=1)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t=1)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{t=1}.
#'
#' CS proximity region is defined with respect to the standard equilateral triangle
#' \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} and edge regions are based on the center of mass \eqn{CM=(1/2,\sqrt{3}/6)}.
#'
#' If \code{p1} and \code{p2} are distinct and either are outside \eqn{T_e}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t=1))} for \code{p1} in \eqn{T_e} that is, returns 1 if \code{p2}
#' is in \eqn{N_{CS}(p1,t=1)}, returns 0 otherwise
#'
#' @seealso \code{\link{IarcCSstd.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#' n<-3
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' IarcCSt1.std.tri(Xp[1,],Xp[2,])
#' IarcCSt1.std.tri(c(.2,.5),Xp[2,])
#' }
#'
#' @export
IarcCSt1.std.tri <- function(p1,p2)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('both arguments must be numeric 2D points')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C)
  if (!in.triangle(p1,Te, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  x<-c(p1[1],p2[1]);y<-c(p1[2],p2[2]);
  if (y[1] <= .5773502694*x[1] && y[1]<= .5773502694-.5773502694*x[1])
  {arc<-IarcCSt1.std.triRAB(p1,p2)}
  else
  {
    if (y[1] > .5773502694-.5773502694*x[1] && x[1] >= .5)
    {arc<-IarcCSt1.std.triRBC(p1,p2)}
    else
    {arc<-IarcCSt1.std.triRAC(p1,p2)}
  }
  arc
} #end of the function
#'

#################################################################

# funsCSEdgeRegs
#'
#' @title Each function is for the presence of an arc from a point in one of the edge regions
#' to another for Central Similarity Proximity Catch Digraphs (CS-PCDs) - standard equilateral triangle case
#'
#' @description
#' Three indicator functions: \code{IarcCSstd.triRAB}, \code{IarcCSstd.triRBC} and \code{IarcCSstd.triRAC}.
#'
#' The function \code{IarcCSstd.triRAB} returns I(\code{p2} is in \eqn{N_{CS}(p1,t)} for \code{p1} in \eqn{RAB} (edge region for edge \eqn{AB},
#' i.e., edge 3) in the standard equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))};
#'
#' \code{IarcCSstd.triRBC} returns I(\code{p2} is in \eqn{N_{CS}(p1,t)} for \code{p1} in \eqn{RBC} (edge region for edge \eqn{BC}, i.e., edge 1) in \eqn{T_e};
#' and
#'
#' \code{IarcCSstd.triRAC} returns I(\code{p2} is in \eqn{N_{CS}(p1,t)} for \code{p1} in \eqn{RAC} (edge region for edge \eqn{AC}, i.e., edge 2) in \eqn{T_e}.
#' That is, each function returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)}, returns 0 otherwise.
#'
#' CS proximity region is defined with respect to \eqn{T_e} whose vertices are also labeled as \eqn{T_e=T(v=1,v=2,v=3)}
#' with expansion parameter \eqn{t>0} and edge regions are based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e}
#'
#' If \code{p1} and \code{p2} are distinct and \code{p1} is outside the corresponding edge region and \code{p2} is outside \eqn{T_e}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their location (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}.
#'
#' @return Each function returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for \code{p1}, that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise
#'
#' @name funsCSEdgeRegs
NULL
#'
#' @seealso \code{\link{IarcCSt1.std.triRAB}}, \code{\link{IarcCSt1.std.triRBC}} and \code{\link{IarcCSt1.std.triRAC}}
#'
#' @rdname funsCSEdgeRegs
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for IarcCSstd.triRAB
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' T3<-rbind(A,B,CM);
#'
#' set.seed(1)
#' Xp<-runif.std.tri(3)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-1
#'
#' IarcCSstd.triRAB(Xp[1,],Xp[2,],t,M)
#' IarcCSstd.triRAB(c(.2,.5),Xp[2,],t,M)
#' }
#'
#' @export
IarcCSstd.triRAB <- function(p1,p2,t,M)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  tri<-rbind(A,B,M)
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  arc<-0
  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  CM<-(A+B+C)/3;

  if (isTRUE(all.equal(M,CM)))
  {
    if (t <=1)
    {
      if (y[2] > y[1]-t*y[1]
          &&
          y[2] < y[1]-1.732050808*x[1]+2*t*y[1]+1.732050808*x[2]
          &&
          y[2] < y[1]+1.732050808*x[1]+2*t*y[1]-1.732050808*x[2])
      {
        arc<-1
      }
    }
    else
    {
      if (y[2] < y[1]-1.732050808*x[1]+2*t*y[1]+1.732050808*x[2]
          &&
          y[2] < y[1]+1.732050808*x[1]+2*t*y[1]-1.732050808*x[2])
      {
        arc<-1
      }
    }
  } else
  {
    reg<-NCStri(p1,Te,t,M,re=3)
    arc<-ifelse(area.polygon(reg)==0,0,sum(in.triangle(p2,reg,boundary=TRUE)$in.tri))
  }
  arc
} #end of the function
#'
#' @rdname funsCSEdgeRegs
#'
#' @examples
#' \dontrun{
#' #Examples for IarcCSstd.triRBC
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' T1<-rbind(B,C,CM);
#'
#' set.seed(1)
#' Xp<-runif.std.tri(3)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-1
#'
#' IarcCSstd.triRBC(Xp[1,],Xp[2,],t,M)
#' IarcCSstd.triRBC(c(.2,.5),Xp[2,],t,M)
#' }
#'
#' @export
IarcCSstd.triRBC <- function(p1,p2,t,M)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  tri<-rbind(B,C,M);
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  arc<-0
  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  CM<-(A+B+C)/3;

  if (isTRUE(all.equal(M,CM)))
  {
    if (t <=1)
    {
      if (y[2] > y[1]-.8660254042*t+.5*t*y[1]+.8660254042*t*x[1]
          &&
          y[2] < y[1]+1.732050808*x[1]-1.732050808*x[2]+1.732050808*t-t*y[1]-1.732050808*t*x[1]
          &&
          y[2] < y[1]-1.732050808*x[1]+1.732050808*x[2]+1.732050808*t-t*y[1]-1.732050808*t*x[1])
      {
        arc<-1
      }
    }
    else
    {
      if (y[2] > y[1]-.8660254042*t+.5*t*y[1]+.8660254042*t*x[1]
          &&
          y[2] < y[1]-1.732050808*x[1]+1.732050808*x[2]+1.732050808*t-t*y[1]-1.732050808*t*x[1])
      {
        arc<-1
      }
    }
  } else
  {
    reg<-NCStri(p1,Te,t,M,re=1)
    arc<-ifelse(area.polygon(reg)==0,0,sum(in.triangle(p2,reg,boundary=TRUE)$in.tri))
  }
  arc
} #end of the function
#'
#' @rdname funsCSEdgeRegs
#'
#' @examples
#' \dontrun{
#' #Examples for IarcCSstd.triRAC
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' T2<-rbind(A,C,CM);
#'
#' set.seed(1)
#' Xp<-runif.std.tri(3)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-1
#'
#' IarcCSstd.triRAC(Xp[1,],Xp[2,],t,M)
#' IarcCSstd.triRAC(c(.2,.5),Xp[2,],t,M)
#' }
#'
#' @export
IarcCSstd.triRAC <- function(p1,p2,t,M)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  tri<-rbind(A,C,M);
  if (!in.triangle(p1,tri, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  arc<-0
  x<-c(p1[1],p2[1]); y<-c(p1[2],p2[2]);
  CM<-(A+B+C)/3;

  if (isTRUE(all.equal(M,CM)))
  {
    if (t <=1)
    {
      if (y[2] > y[1]+.5*t*y[1]-.8660254042*t*x[1]
          &&
          y[2] < y[1]-1.732050808*x[1]+1.732050808*x[2]-t*y[1]+1.732050808*t*x[1]
          &&
          y[2] < y[1]+1.732050808*x[1]-1.732050808*x[2]-t*y[1]+1.732050808*t*x[1])
      {
        arc<-1
      }
    }
    else
    {
      if (y[2] > y[1]+.5*t*y[1]-.8660254042*t*x[1]
          &&
          y[2] < y[1]+1.732050808*x[1]-1.732050808*x[2]-t*y[1]+1.732050808*t*x[1])
      {
        arc<-1
      }
    }
  } else
  {
    reg<-NCStri(p1,Te,t,M,re=2)
    arc<-ifelse(area.polygon(reg)==0,0,sum(in.triangle(p2,reg,boundary=TRUE)$in.tri))
  }
  arc
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for Central Similarity Proximity Catch
#' Digraphs (CS-PCDs) - standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with expansion parameter \eqn{t>0}.
#' This function is equivalent to \code{IarcCSstd.tri}, except that it computes the indicator using the functions
#' \code{IarcCSstd.triRAB}, \code{IarcCSstd.triRBC} and \code{IarcCSstd.triRAC} which are edge-region specific indicator functions.
#' For example,
#' \code{IarcCSstd.triRAB} computes \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2} when \code{p1}
#' resides in the edge region of edge \eqn{AB}.
#'
#' CS proximity region is defined with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e}.
#' \code{re} is the index of the edge region \code{p1} resides, with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \eqn{T_e}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param re The index of the edge region in \eqn{T_e} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for \code{p1},
#' that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)}, returns 0 otherwise
#'
#' @seealso \code{\link{IarcCStri}} and \code{\link{IarcPEstd.tri}}
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
#' n<-3
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-1
#' IarcCSedge.reg.std.tri(Xp[1,],Xp[2,],t,M)
#' IarcCSstd.tri(Xp[1,],Xp[2,],t,M)
#'
#' #or try
#' re<-rel.edge.std.triCM(Xp[1,])$re
#' IarcCSedge.reg.std.tri(Xp[1,],Xp[2,],t,M,re=re)
#' }
#'
#' @export IarcCSedge.reg.std.tri
IarcCSedge.reg.std.tri <- function(p1,p2,t,M=c(1,1,1),re=NULL)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,Te, boundary = TRUE)$in.tri || !in.triangle(p2,Te, boundary = TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(re))
  {re<-rel.edge.tri(p1,Te,M)$re #edge region for p1
  } else
  {  if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')}}

  x<-c(p1[1],p2[1]);y<-c(p1[2],p2[2]);
  if (re==3)
  {arc<-IarcCSstd.triRAB(p1,p2,t,M)
  } else
  {
    if (re==1)
    {arc<-IarcCSstd.triRBC(p1,p2,t,M)}
    else
    {arc<-IarcCSstd.triRAC(p1,p2,t,M)}
  }
  arc
} #end of the function
#'

#################################################################

#' @title Incidence matrix for Central Similarity Proximity Catch Digraphs (CS-PCDs) - standard
#' equilateral triangle case
#'
#' @description Returns the incidence matrix for the CS-PCD whose vertices are the given 2D numerical data set, \code{Xp},
#' in the standard equilateral triangle \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#'
#' CS proximity region is defined with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates.
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return Incidence matrix for the CS-PCD with vertices being 2D data set, \code{Xp} and CS proximity
#' regions are defined in the standard equilateral triangle \eqn{T_e} with \code{M}-edge regions.
#'
#' @seealso \code{\link{inci.matCStri}}, \code{\link{inci.matCS}} and \code{\link{inci.matPEstd.tri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' inc.mat<-inci.matCSstd.tri(Xp,t=1.25,M)
#' inc.mat
#' sum(inc.mat)-n
#' num.arcsCSstd.tri(Xp,t=1.25)
#'
#' dom.num.greedy(inc.mat) #try also dom.num.exact(inc.mat)  #might take a long time for large n
#' Idom.num.up.bnd(inc.mat,1)
#' }
#'
#' @export inci.matCSstd.tri
inci.matCSstd.tri <- function(Xp,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp)

  inc.mat<-matrix(0, nrow=n, ncol=n)
  for (i in 1:n)
  {p1<-Xp[i,]
  for (j in ((1:n)) )
  {p2<-Xp[j,]
  inc.mat[i,j]<-IarcCSstd.tri(p1,p2,t,M)
  }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title Number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' and quantities related to the triangle - standard equilateral triangle case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}.
#' It also provides number of vertices
#' (i.e., number of data points inside the standard equilateral triangle \eqn{T_e})
#' and indices of the data points that reside in \eqn{T_e}.
#'
#' CS proximity region \eqn{N_{CS}(x,t)} is defined with respect to the standard
#' equilateral triangle \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with expansion parameter \eqn{t>0}
#' and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e}.
#' For the number of arcs, loops are not allowed so
#' arcs are only possible for points inside \eqn{T_e} for this function.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the digraph.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates.
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and quantities related to the standard equilateral triangle}
#' \item{num.arcs}{Number of arcs of the CS-PCD}
#' \item{num.in.tri}{Number of \code{Xp} points in the standard equilateral triangle, \eqn{T_e}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points that reside in \eqn{T_e}}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support triangle \eqn{T_e}.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsCStri}}, \code{\link{num.arcsCS}}, and \code{\link{num.arcsPEstd.tri}},
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' Narcs = num.arcsCSstd.tri(Xp,t=.5,M)
#' Narcs
#' summary(Narcs)
#' par(pty="s")
#' plot(Narcs,asp=1)
#' }
#'
#' @export num.arcsCSstd.tri
num.arcsCSstd.tri <- function(Xp,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp)
  arcs<-0
  ind.in.tri = NULL
  if (n<=0)
  {
    arcs<-0
  } else
  {
  for (i in 1:n)
  {p1<-as.numeric(Xp[i,])
  if (!in.triangle(p1,Te,boundary = TRUE)$in.tri) {
    arcs<-arcs+0
  } else
  {
    ind.in.tri = c(ind.in.tri,i)

    for (j in ((1:n)[-i]) ) #to avoid loops
    {p2<-as.numeric(Xp[j,])
    if (!in.triangle(p2,Te,boundary = TRUE)$in.tri)
    {arcs<-arcs+0
    } else
    {
      arcs<-arcs + IarcCSstd.tri(p1,p2,t,M)
    }
    }
  }
  }
}

  NinTri = length(ind.in.tri)
  desc<-"Number of Arcs of the CS-PCD and the Related Quantities with vertices Xp in the Standard Equilateral Triangle"

  res<-list(desc=desc, #description of the output
            num.arcs=arcs, #number of arcs for the AS-PCD
            num.in.tri=NinTri, # number of Xp points in CH of Yp points
            ind.in.tri=ind.in.tri, #indices of data points inside the triangle
            tess.points=Te, #tessellation points
            vertices=Xp #vertices of the digraph
  )

  class(res)<-"NumArcs"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' and quantities related to the triangle - one triangle case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' CS proximity region \eqn{N_{CS}(x,t)} is defined with respect to the triangle, \code{tri}  with expansion parameter \eqn{t>0}
#' and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' For the number of arcs, loops are not allowed so
#' arcs are only possible for points inside \code{tri} for this function.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}; default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and quantities related to the triangle}
#' \item{num.arcs}{Number of arcs of the CS-PCD}
#' \item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points that reside in the triangle}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support triangle.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsCSstd.tri}}, \code{\link{num.arcsCS}}, \code{\link{num.arcsPEtri}},
#' and \code{\link{num.arcsAStri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' Narcs = num.arcsCStri(Xp,Tr,t=.5,M)
#' Narcs
#' summary(Narcs)
#' plot(Narcs)
#' }
#'
#' @export num.arcsCStri
num.arcsCStri <- function(Xp,tri,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or
        3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp)
  arcs<-0
  ind.in.tri = NULL
  if (n<=0)
  {
    arcs<-0
  } else
  {
    for (i in 1:n)
    { Xpi<-as.numeric(Xp[i,])
    if (in.triangle(Xpi,tri,boundary=TRUE)$in.tri)
    {   edgei<-rel.edges.tri(Xpi,tri,M)$re
    ind.in.tri = c(ind.in.tri,i)
    for (j in (1:n)[-i])  #to avoid loops
    { Xpj<-as.numeric(Xp[j,])
    arcs<-arcs+IarcCStri(Xpi,Xpj,tri,t,M,re=edgei)
    }
    }
    }
  }
  NinTri = length(ind.in.tri)
  desc<-"Number of Arcs of the CS-PCD with vertices Xp and Quantities Related to the Support Triangle"

  res<-list(desc=desc, #description of the output
            num.arcs=arcs, #number of arcs for the AS-PCD
            num.in.tri=NinTri, # number of Xp points in CH of Yp points
            ind.in.tri=ind.in.tri, #indices of data points inside the triangle
            tess.points=tri, #tessellation points
            vertices=Xp #vertices of the digraph
  )

  class(res)<-"NumArcs"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Number of arcs of Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' and related quantities of the induced subdigraphs for points in the Delaunay triangles -
#' multiple triangle case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs and various other quantities related to the Delaunay triangles
#' for Central Similarity Proximity Catch Digraph
#' (CS-PCD) whose vertices are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' CS proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points with expansion parameter \eqn{t>0} and edge regions in each triangle
#' is based on the center \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of each
#' Delaunay triangle or based on circumcenter of each Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
#' Each Delaunay triangle is first converted to an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles whose union constitutes the
#' convex hull of \code{Yp} points). For the number of arcs, loops are not allowed so arcs are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}) for more on CS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates which serves as a center in the interior of each Delaunay
#' triangle, default for \eqn{M=(1,1,1)} which is the center of mass of each triangle.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and related quantities for the induced subdigraphs in the Delaunay triangles}
#' \item{num.arcs}{Total number of arcs in all triangles, i.e., the number of arcs for the entire PE-PCD}
#' \item{num.in.conhull}{Number of \code{Xp} points in the convex hull of \code{Yp} points}
#' \item{num.in.tris}{The vector of number of \code{Xp} points in the Delaunay triangles based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the areas of Delaunay triangles based on \code{Yp} points}
#' \item{tri.num.arcs}{The \code{vector} of the number of arcs of the component of the PE-PCD in the
#' Delaunay triangles based on \code{Yp} points}
#' \item{del.tri.ind}{A matrix of indices of vertices of the Delaunay triangles based on \code{Yp} points,
#' each column corresponds to the vector of indices of the vertices of one triangle.}
#' \item{data.tri.ind}{A \code{vector} of indices of vertices of the Delaunay triangles in which data points reside,
#' i.e., column number of \code{del.tri.ind} for each \code{Xp} point.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the Delaunay triangulation based on \code{Yp} points.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsCStri}}, \code{\link{num.arcsCSstd.tri}}, \code{\link{num.arcsPE}},
#' and \code{\link{num.arcsAS}}
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
#' Xp<-cbind(runif(nx),runif(nx))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' Narcs = num.arcsCS(Xp,Yp,t=1,M)
#' Narcs
#' summary(Narcs)
#' plot(Narcs)
#' }
#'
#' @export num.arcsCS
num.arcsCS <- function(Xp,Yp,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M,3) || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates')}

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
  { tri<-as.basic.tri(Yp)$tri
  NumArcs = num.arcsCStri(Xdt,tri,t,M)
  NinTri<-NumArcs$num.in.tri #number of points in the triangle

  if (NinTri==0)
  {Tot.Arcs<-0;
  ni.vec<-arcs<-rep(0,ndt)
  data.tri.ind = NA
  } else
  {
    Xdt<-matrix(Xp[NumArcs$ind.in.tri,],ncol=2)
    tri<-as.basic.tri(Yp)$tri #convert the triangle Yp into an nonscaled basic triangle, see as.basic.tri help page
    Wvec<-area.polygon(tri)
    Tot.Arcs<- NumArcs$num.arcs #number of arcs in the triangle Yp
    ni.vec = NumArcs$num.in.tri
    Tri.Ind = NumArcs$ind.in.tri
    data.tri.ind = rep(NA,nx)
    data.tri.ind[Tri.Ind] =1
    arcs = NumArcs$num.arcs
  }

  desc<-"Number of Arcs of the CS-PCD with vertices Xp and Related Quantities for the Induced Subdigraph for the Points in the Delaunay Triangles"
  res<-list(desc=desc, #description of the output
            num.arcs=Tot.Arcs,
            tri.num.arcs=arcs,
            num.in.conv.hull=NinTri,
            num.in.tris=ni.vec,
            weight.vec=Wvec,
            del.tri.ind=t(Ytri),
            data.tri.ind=data.tri.ind,
            tess.points=Yp, #tessellation points
            vertices=Xp #vertices of the digraph
  )
  } else
  {
    if (NinCH==0)
    {Tot.Arcs<-0;
    ni.vec<-arcs<-rep(0,ndt)
    data.tri.ind =NULL
    } else
    {
      Tri.Ind<-indices.delaunay.tri(Xp,Yp,Ytrimesh) #indices of triangles in which the points in the data fall
      ind.in.CH = which(!is.na(Tri.Ind))
      #calculation of the total number of arcs
      ni.vec<-arcs<-vector()
      data.tri.ind = rep(NA,nx)
      for (i in 1:ndt)
      {
        dt.ind=which(Tri.Ind==i) #which indices of data points residing in ith Delaunay triangle
        Xpi<-Xp[dt.ind,] #points in ith Delaunay triangle
        data.tri.ind[dt.ind] =i #assigning the index of the Delaunay triangle that contains the data point
        ifelse(ndt==1,Tri<-Yp[Ytri,],Tri<-Yp[Ytri[i,],])  #vertices of ith triangle
        tri<-as.basic.tri(Tri)$tri #convert the triangle Tri into an nonscaled basic triangle, see as.basic.tri help page
        ni.vec<-c(ni.vec,length(Xpi)/2)  #number of points in ith Delaunay triangle

        num.arcs<-num.arcsCStri(Xpi,tri,t,M)$num.arcs  #number of arcs in ith triangle
        arcs<-c(arcs,num.arcs)  #number of arcs in all triangles as a vector

      }

      Tot.Arcs<-sum(arcs)  #the total number of arcs in all triangles
    }

    desc<-"Number of Arcs of the CS-PCD with vertices Xp and Related Quantities for the Induced Subdigraphs for the Points in the Delaunay Triangles"

    res<-list(desc=desc, #description of the output
              num.arcs=Tot.Arcs, #number of arcs for the entire PCD
              tri.num.arcs=arcs, #vector of number of arcs for the Delaunay triangles
              num.in.conv.hull=NinCH, #number of Xp points in CH of Yp points
              ind.in.conv.hull= ind.in.CH, #indices of Xp points in CH of Yp points
              num.in.tris=ni.vec, # vector of number of Xp points in the Delaunay triangles
              weight.vec=Wvec, #areas of Delaunay triangles
              del.tri.ind=t(Ytri), # indices of the vertices of the Delaunay triangles, i.e., each column corresponds to the indices of the vertices of one Delaunay triangle
              data.tri.ind=data.tri.ind, #indices of the Delaunay triangles in which data points reside, i.e., column number of del.tri.ind for each Xp point
              tess.points=Yp, #tessellation points
              vertices=Xp #vertices of the digraph
    )
  }

  class(res)<-"NumArcs"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title A test of segregation/association based on arc density of Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 2D data
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points in the convex hull of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster away from \code{Yp} points) and association (where \code{Xp} points cluster around
#' \code{Yp} points) based on the normal approximation of the arc density of the CS-PCD for uniform 2D data
#' in the convex hull of \code{Yp} points.
#'
#' The function yields the test statistic, \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval, estimate and null value for the parameter of interest (which is the arc density),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points in the convex hull of \code{Yp} points, arc density
#' of CS-PCD whose vertices are \code{Xp} points equals to its expected value under the uniform distribution and
#' \code{alternative} could be two-sided, or left-sided (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated around the centers of the triangles, or segregation).
#'
#' CS proximity region is constructed with the expansion parameter \eqn{t>0} and \eqn{CM}-edge regions
#' (i.e., the test is not available for a general center \eqn{M} at this version of the function).
#'
#' **Caveat:** This test is currently a conditional test, where \code{Xp} points are assumed to be random, while \code{Yp} points are
#' assumed to be fixed (i.e., the test is conditional on \code{Yp} points).
#' Furthermore, the test is a large sample test when \code{Xp} points are substantially larger than \code{Yp} points,
#' say at least 5 times more.
#' This test is more appropriate when supports of \code{Xp} and \code{Yp} has a substantial overlap.
#' Currently, the \code{Xp} points outside the convex hull of \code{Yp} points are handled with a convex hull correction factor
#' (see the description below and the function code.)
#' However, in the special case of no \code{Xp} points in the convex hull of \code{Yp} points, arc density is taken to be 1,
#' as this is clearly a case of segregation. Removing the conditioning and extending it to the case of non-concurring supports is
#' an ongoing line of research of the author of the package.
#'
#' \code{ch.cor} is for convex hull correction (default is \code{"no convex hull correction"}, i.e., \code{ch.cor=FALSE})
#' which is recommended when both \code{Xp} and \code{Yp} have the same rectangular support.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param ch.cor A logical argument for convex hull correction, default \code{ch.cor=FALSE},
#' recommended when both \code{Xp} and \code{Yp} have the same rectangular support.
#' @param alternative Type of the alternative hypothesis in the test, one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval, default is \code{0.95}, for the arc density of CS-PCD based on
#' the 2D data set \code{Xp}.
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
#' @seealso \code{\link{PEarc.dens.test}} and \code{\link{CSarc.dens.test1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-5;  #try also nx<-40; ny<-10 or nx<-1000; ny<-10;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx),runif(nx))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' plotDelaunay.tri(Xp,Yp,xlab="",ylab = "")
#'
#' CSarc.dens.test(Xp,Yp,t=.5)
#' CSarc.dens.test(Xp,Yp,t=.5,ch=TRUE)
#' #try also t=1.0 and 1.5 above
#' }
#'
#' @export CSarc.dens.test
CSarc.dens.test <- function(Xp,Yp,t,ch.cor=FALSE,
                            alternative = c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one \"greater\", \"less\", \"two.sided\"")

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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  Arcs<-num.arcsCS(Xp,Yp,t,M=c(1,1,1))  #use the default, i.e., CM for the center M
  NinCH<-Arcs$num.in.con

  num.arcs<-Arcs$num.arcs #total number of arcs in the PE-PCD
  num.arcs.tris = Arcs$tri.num.arcs #vector of number of arcs in the Delaunay triangles
  num.dat.tris = Arcs$num.in.tris #vector of number of data points in the Delaunay triangles
  Wvec<-Arcs$w
  LW<-Wvec/sum(Wvec)

  tri.ind = Arcs$data.tri.ind
  ind.triCH =  t(Arcs$del.tri)

  ind.Xp1 = which(num.dat.tris==1)
  if (length(ind.Xp1)>0)
  {
    for (i in ind.Xp1)
    {
      Xpi = Xp[which(tri.ind==i),]
      tri =  Yp[ind.triCH[i,],]
      npe = NCStri(Xpi,tri,t)
      num.arcs = num.arcs+area.polygon(npe)/Wvec[i]
    }
  }

  asy.mean0<-muCS2D(t)  #asy mean value for the t value
  asy.mean<-asy.mean0*sum(LW^2)

  asy.var0<-asyvarCS2D(t)  #asy variance value for the t value
  asy.var<-asy.var0*sum(LW^3)+4*asy.mean0^2*(sum(LW^3)-(sum(LW^2))^2)

  n<-nrow(Xp)  #number of X points
  if (NinCH == 0)
  {warning('There is no Xp point in the convex hull of Yp points to compute arc density,
           but as this is clearly a segregation pattern, arc density is taken to be 1!')
    arc.dens=1
    TS0<-sqrt(n)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat
  } else
  {  arc.dens<-num.arcs/(NinCH*(NinCH-1))
  TS0<-sqrt(NinCH)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat}  #arc density
  }
  estimate1<-arc.dens
  estimate2<-asy.mean

  method <-c("Large Sample z-Test Based on Arc Density of CS-PCD for Testing Uniformity of 2D Data ---")
  if (ch.cor==FALSE)
  {
    TS<-TS0
    method <-c(method, "\n without Convex Hull Correction")
  }
  else
  {
    m<-nrow(Yp)  #number of Y points
    NoutCH<-n-NinCH #number of points outside of the convex hull

    prop.out <- NoutCH/n #observed proportion of points outside convex hull
    exp.prop.out <- 1.7932/m+1.2229/sqrt(m) #expected proportion of points outside convex hull

    TS<-TS0+abs(TS0)*sign(prop.out-exp.prop.out)*(prop.out-exp.prop.out)^2
    method <-c(method, "\n with Convex Hull Correction")
  }

  names(estimate1) <-c("arc density")
  null.dens<-asy.mean
  names(null.dens) <-"(expected) arc density"
  names(TS) <-"standardized arc density (i.e., Z)"

  if (alternative == "less") {
    pval <-pnorm(TS)
    cint <-arc.dens+c(-Inf, qnorm(conf.level))*sqrt(asy.var/NinCH)
  }
  else if (alternative == "greater") {
    pval <-pnorm(TS, lower.tail = FALSE)
    cint <-arc.dens+c(-qnorm(conf.level),Inf)*sqrt(asy.var/NinCH)
  }
  else {
    pval <-2 * pnorm(-abs(TS))
    alpha <-1 - conf.level
    cint <-qnorm(1 - alpha/2)
    cint <-arc.dens+c(-cint, cint)*sqrt(asy.var/NinCH)
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

#' @title The arcs of Central Similarity Proximity Catch Digraphs (CS-PCD) for 2D data - one triangle case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) for data set \code{Xp} as the vertices
#' of CS-PCD and related parameters and the quantities of the digraph.
#'
#' CS proximity regions are constructed with respect to the triangle \code{tri} with expansion
#' parameter \eqn{t>0}, i.e., arcs may exist for points only inside \code{tri}.
#' It also provides various descriptions and quantities about the arcs of the CS-PCD
#' such as number of arcs, arc density, etc.
#'
#' Edge regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of
#' the triangle \code{tri}; default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri} or the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, the center \code{M} used to
#' construct the edge regions and the expansion parameter \code{t}.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support triangle.}
#' \item{tess.name}{Name of data set used in tessellation (i.e., vertices of the triangle)}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of CS-PCD for 2D data set \code{Xp} as vertices of the digraph}
#' \item{E}{Heads (or arrow ends) of the arcs of CS-PCD for 2D data set \code{Xp} as vertices of the digraph}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of triangles, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsCS}}, \code{\link{arcsAStri}} and \code{\link{arcsPEtri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' t<-1.5  #try also t<-2
#'
#' Arcs<-arcsCStri(Xp,Tr,t,M)
#' #or try with the default center Arcs<-arcsCStri(Xp,Tr,t); M= (Arcs$param)$c
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' #can add edge regions
#' L<-rbind(M,M,M); R<-Tr
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' #now we can add the vertex names and annotation
#' txt<-rbind(Tr,M)
#' xc<-txt[,1]+c(-.02,.03,.02,.03)
#' yc<-txt[,2]+c(.02,.02,.03,.06)
#' txt.str<-c("A","B","C","M")
#' text(xc,yc,txt.str)
#' }
#'
#' @export arcsCStri
arcsCStri <- function(Xp,tri,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp)

  in.tri<-rep(0,n)
  for (i in 1:n)
    {in.tri[i]<-in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri } #indices the Xp points inside the triangle

  Xtri<-Xp[in.tri==1,] #the Xp points inside the triangle
  n2<-length(Xtri)/2

  #the arcs of CS-PCDs
  S<-E<-NULL #S is for source and E is for end points for the arcs
  if (n2>1)
  {
    for (j in 1:n2)
    {
      p1<-Xtri[j,]; RE1<-rel.edge.tri(p1,tri,M)$re;
      for (k in (1:n2)[-j])  #to avoid loops
      {
        p2<-Xtri[k,];
       # print(IarcCStri(p1,p2,tri,t,M,RE1))
        if (IarcCStri(p1,p2,tri,t,M,RE1)==1)
        {
          S <-rbind(S,Xtri[j,]); E <-rbind(E,Xtri[k,]);
        }
      }
    }
  }

  cname <-"M"
  param<-list(M,t)
  names(param)<-c("center","expansion parameter")

  Mr<-round(M,2)
  main.txt<-paste("Arcs of CS-PCD\n with t = ",t," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
  typ<-paste("Central Similarity Proximity Catch Digraph (CS-PCD) for 2D Points in the Triangle with Expansion Parameter t = ",t," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")

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

#' @title Incidence matrix for Central Similarity Proximity Catch Digraphs (CS-PCDs) - one triangle case
#'
#' @description Returns the incidence matrix for the CS-PCD whose vertices are the given 2D numerical data set, \code{Xp},
#' in the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)}.
#'
#' CS proximity regions are constructed with respect to triangle \code{tri}
#' with expansion parameter \eqn{t>0} and edge regions are based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' @return Incidence matrix for the CS-PCD with vertices being 2D data set, \code{Xp},
#' in the triangle \code{tri} with edge regions based on center \code{M}
#'
#' @seealso \code{\link{inci.matCS}}, \code{\link{inci.matPEtri}}, and  \code{\link{inci.matAStri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' IM<-inci.matCStri(Xp,Tr,t=1.25,M)
#' IM
#' dom.num.greedy(IM) #try also dom.num.exact(IM)
#' Idom.num.up.bnd(IM,3)
#'
#' inci.matCStri(Xp,Tr,t=1.5,M)
#' }
#'
#' @export inci.matCStri
inci.matCStri <- function(Xp,tri,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (isTRUE(all.equal(M,circumcenter.tri(tri)))==FALSE & in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  n<-nrow(Xp)

  inc.mat<-matrix(0, nrow=n, ncol=n)
  if (n>1)
  {
    for (i in 1:n)
    {p1<-Xp[i,]
    RE<-rel.edge.tri(p1,tri,M)$re
    for (j in 1:n )
    {p2<-Xp[j,]
    inc.mat[i,j]<-IarcCStri(p1,p2,tri,t,M,re=RE)
    }
    }
  }
  diag(inc.mat)<-1
  inc.mat
} #end of the function
#'

#################################################################

#' @title The plot of the arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for a
#' 2D data set - one triangle case
#'
#' @description Plots the arcs of CS-PCD whose vertices are the data points, \code{Xp} and the triangle \code{tri}. CS proximity regions
#' are constructed with respect to the triangle \code{tri} with expansion parameter \eqn{t>0}, i.e., arcs may exist only
#' for \code{Xp} points inside the triangle \code{tri}.
#'
#' Edge regions are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}; default
#' is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' @param asp A \code{numeric} value, giving the aspect ratio \eqn{y/x} (default is \code{NA}), see the official help page for \code{asp} by
#' typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param edge.reg A logical argument to add edge regions to the plot, default is \code{edge.reg=FALSE}.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of the CS-PCD whose vertices are the points in data set \code{Xp} and the triangle \code{tri}
#'
#' @seealso \code{\link{plotCSarcs}}, \code{\link{plotPEarcs.tri}} and \code{\link{plotASarcs.tri}}
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
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' t<-1.5  #try also t<-2
#'
#' plotCSarcs.tri(Xp,Tr,t,M,main="Arcs of CS-PCD with t=1.5",xlab="",ylab="",edge.reg = TRUE)
#'
#' # or try the default center
#' #plotCSarcs.tri(Xp,Tr,t,main="Arcs of CS-PCD with t=1.5",xlab="",ylab="",edge.reg = TRUE);
#' #M=(arcsCStri(Xp,Tr,r)$param)$c #the part "M=(arcsPEtri(Xp,Tr,r)$param)$cent" is optional,
#' #for the below annotation of the plot
#'
#' #can add vertex labels and text to the figure (with edge regions)
#' txt<-rbind(Tr,M)
#' xc<-txt[,1]+c(-.02,.02,.02,.03)
#' yc<-txt[,2]+c(.02,.02,.02,.03)
#' txt.str<-c("A","B","C","M")
#' text(xc,yc,txt.str)
#' }
#'
#' @export plotCSarcs.tri
plotCSarcs.tri <- function(Xp,tri,t,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,
                           xlim=NULL,ylim=NULL,edge.reg=FALSE,...)
{
  arcsCS<-arcsCStri(Xp,tri,t,M)
  S<-arcsCS$S
  E<-arcsCS$E
  cent = (arcsCS$param)$c

  Xp<-matrix(Xp,ncol=2)
  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2])}

  Mr=round(cent,2)
  if (is.null(main))
  {Mvec= paste(Mr, collapse=",")
  main=paste("Arcs of CS-PCD\n with t = ",t," and M = (",Mvec,")",sep="")
  }

  if (edge.reg)
  {main=c(main,"\n (edge regions added)")}

  plot(Xp,main=main,asp=asp, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3,...)
  polygon(tri)
  if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}

  if (edge.reg){
    L<-rbind(cent,cent,cent); R<-tri
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }
} #end of the function
#'

#################################################################

#' @title The plot of the Central Similarity (CS) Proximity Regions for a 2D data set - one triangle case
#'
#' @description Plots the points in and outside of the triangle \code{tri} and also the CS proximity regions which are also
#' triangular for points inside the triangle \code{tri} with edge regions are based on the center of mass CM.
#'
#' CS proximity regions are defined with respect to the triangle \code{tri}
#' with expansion parameter \eqn{t>0}, so CS proximity regions are defined only for points inside the
#' triangle \code{tri}.
#'
#' Edge regions are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param Xp A set of 2D points for which CS proximity regions are constructed.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' @param asp A \code{numeric} value, giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param edge.reg A logical argument to add edge regions to the plot, default is \code{edge.reg=FALSE}.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the CS proximity regions for points inside the triangle \code{tri}
#' (and just the points outside \code{tri})
#'
#' @seealso \code{\link{plotCSregs}}, \code{\link{plotASregs.tri}} and \code{\link{plotPEregs.tri}},
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' t<-.5  #try also t<-2
#'
#' plotCSregs.tri(Xp0,Tr,t,M,main="Proximity Regions for CS-PCD", xlab="",ylab="")
#'
#' Xp = Xp0[1,]
#' plotCSregs.tri(Xp,Tr,t,M,main="CS Proximity Regions with t=.5", xlab="",ylab="",edge.reg=TRUE)
#'
#' # or try the default center
#' plotCSregs.tri(Xp,Tr,t,main="CS Proximity Regions with t=.5", xlab="",ylab="",edge.reg=TRUE);
#' #M=(arcsCStri(Xp,Tr,r)$param)$c #the part "M=(arcsPEtri(Xp,Tr,r)$param)$cent" is optional,
#' #for the below annotation of the plot
#'
#' #can add vertex labels and text to the figure (with edge regions)
#' txt<-rbind(Tr,M)
#' xc<-txt[,1]+c(-.02,.02,.02,.02)
#' yc<-txt[,2]+c(.02,.02,.02,.03)
#' txt.str<-c("A","B","C","M")
#' text(xc,yc,txt.str)
#' }
#'
#' @export plotCSregs.tri
plotCSregs.tri <- function(Xp,tri,t,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,edge.reg=FALSE,...)
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp)

  in.tri<-rep(0,n)
  for (i in 1:n)
    in.tri[i]<-in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri #indices of the Xp points inside the triangle

  Xtri<-matrix(Xp[in.tri==1,],ncol=2)  #the Xp points inside the triangle
  nt<-length(Xtri)/2 #number of Xp points inside the triangle

  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2])}

  xr<-xlim[2]-xlim[1]
  yr<-ylim[2]-ylim[1]

  Mr=round(M,2)
  if (is.null(main))
  {Mvec= paste(Mr, collapse=",")
  main=paste("CS Proximity Regions\n with t = ",t," and M = (",Mvec,")",sep="")
  }

  if (edge.reg)
  {main=c(main,"\n (edge regions added)")}

  plot(Xp,main=main, asp=asp, xlab=xlab, ylab=ylab,xlim=xlim+xr*c(-.05,.05),
       ylim=ylim+yr*c(-.05,.05),pch=".",cex=3,...)
  polygon(tri,lty=2)
  if (nt>=1)
  {
    for (i in 1:nt)
    {
      P1<-Xtri[i,]
      RE<-rel.edge.tri(P1,tri,M)$re

      pr<-NCStri(P1,tri,t,M,re=RE)
      polygon(pr,border="blue")
    }
  }

  if (edge.reg){
    L<-rbind(M,M,M); R<-tri
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }

} #end of the function
#'

#################################################################

#' @title The arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for 2D data - multiple triangle case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) of Central Similarity Proximity Catch Digraph
#' (CS-PCD) whose vertices are the data points in \code{Xp} in the multiple triangle case
#' and related parameters and the quantities of the digraph.
#'
#' CS proximity regions are
#' defined with respect to the Delaunay triangles based on \code{Yp} points with expansion parameter \eqn{t>0} and
#' edge regions in each triangle are based on the center \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of each Delaunay triangle (default for \eqn{M=(1,1,1)} which is the center of mass of
#' the triangle). Each Delaunay triangle is first converted to an (nonscaled) basic triangle so that
#' \code{M} will be the same type of center for each Delaunay triangle (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles whose union constitutes the
#' convex hull of \code{Yp} points). For the number of arcs, loops are not allowed so arcs are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}) for more on CS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates which serves as a center in the interior of each Delaunay
#' triangle, default for \eqn{M=(1,1,1)} which is the center of mass of each triangle.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph, here, it is the center used to construct the edge regions.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is Delaunay triangulation based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation, it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of CS-PCD for 2D data set \code{Xp} as vertices of the digraph}
#' \item{E}{Heads (or arrow ends) of the arcs of CS-PCD for 2D data set \code{Xp} as vertices of the digraph}
#' \item{mtitle}{Text for \code{"main"} title in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices, number of partition points,
#' number of triangles, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsCStri}}, \code{\link{arcsAS}} and \code{\link{arcsPE}}
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
#' tau<-1.5  #try also tau<-2
#'
#' Arcs<-arcsCS(Xp,Yp,tau,M)
#' #or use the default center Arcs<-arcsCS(Xp,Yp,tau)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#' }
#'
#' @export arcsCS
arcsCS <- function(Xp,Yp,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (nrow(Yp)==3)
  {
    res<-arcsCStri(Xp,Yp,t,M)
  } else
  { if (!is.point(M,3) || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)  #the Xp points inside the convex hull of Yp

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)
    nx2<-nrow(Xch)

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
          Yi.tri<-as.basic.tri(Yi.Tri)$tri #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.basic.tri help page

          nl<-nrow(Xl)
          re.ind<-rel.edges.tri(Xl,Yi.tri,M)$re
          for (j in 1:nl)
          {RE<-re.ind[j]
          for (k in (1:nl)[-j])  #to avoid loops
            if (IarcCStri(Xl[j,],Xl[k,],Yi.tri,t,M,re=RE)==1 )
            {
              S <-rbind(S,Xl[j,]); E <-rbind(E,Xl[k,]);
            }
          }
        }
      }
    }

    cname <-"M"

    param<-list(M,t)
    names(param)<-c("center","expansion parameter")

    Mvec= paste(M, collapse=",")
    typ<-paste("Central Similarity Proximity Catch Digraph (CS-PCD) for 2D Points in the Multiple Triangles with Expansion Parameter t = ",t," and Center ", cname," = (",Mvec,")",sep="")

    main.txt=paste("Arcs of CS-PCD\n with t = ",t," and Center ", cname," = (",Mvec,")",sep="")

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

#' @title Incidence matrix for Central Similarity Proximity Catch Digraphs (CS-PCDs) - multiple triangle case
#'
#' @description Returns the incidence matrix of Central Similarity Proximity Catch Digraph (CS-PCD) whose vertices are the
#' data points in \code{Xp} in the multiple triangle case.
#'
#' CS proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points with expansion parameter \eqn{t>0} and edge regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of each Delaunay
#' triangle (default for \eqn{M=(1,1,1)} which is the center of mass of the triangle).
#' Each Delaunay triangle is first converted to an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles whose union constitutes the
#' convex hull of \code{Yp} points). For the incidence matrix loops are allowed,
#' so the diagonal entries are all equal to 1.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}) for more on CS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates which serves as a center in the interior of each Delaunay
#' triangle, default for \eqn{M=(1,1,1)} which is the center of mass of each triangle.
#'
#' @return Incidence matrix for the CS-PCD with vertices being 2D data set, \code{Xp}.
#' CS proximity regions are constructed with respect to the Delaunay triangles and \code{M}-edge regions.
#'
#' @seealso \code{\link{inci.matCStri}}, \code{\link{inci.matCSstd.tri}}, \code{\link{inci.matAS}},
#' and \code{\link{inci.matPE}}
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
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' t<-1.5  #try also t<-2
#'
#' IM<-inci.matCS(Xp,Yp,t,M)
#' IM
#' dom.num.greedy(IM) #try also dom.num.exact(IM)  #takes a very long time for large nx, try smaller nx
#' Idom.num.up.bnd(IM,3)  #takes a very long time for large nx, try smaller nx
#' }
#'
#' @export inci.matCS
inci.matCS <- function(Xp,Yp,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (nrow(Yp)==3)
  {
    inc.mat<-inci.matCStri(Xp,Yp,t,M)
  } else
  { if (!is.point(M,3) || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates')}

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
        Yi.Tri<-Yp[DTr[i.tr[i],],] #vertices of the ith triangle
        Yi.tri<-as.basic.tri(Yi.Tri)$tri #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.basic.tri help page

        edge<-rel.edge.tri(p1,Yi.tri,M)$re
        for (j in 1:nx )
        {p2<-Xp[j,]
        inc.mat[i,j]<-IarcCStri(p1,p2,Yi.tri,t,M,re=edge)
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

#' @title The plot of the arcs of Central Similarity Proximity Catch Digraph (CS-PCD) for a
#' 2D data set - multiple triangle case
#'
#' @description Plots the arcs of Central Similarity Proximity Catch Digraph (CS-PCD) whose vertices are the data
#' points in \code{Xp} in the multiple triangle case and the Delaunay triangles based on \code{Yp} points.
#'
#' CS proximity regions are defined with respect to the Delaunay triangles based on \code{Yp} points with
#' expansion parameter \eqn{t>0} and edge regions in each triangle are based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle). Each Delaunay triangle is first converted to an (nonscaled)
#' basic triangle so that \code{M} will be the same type of center for each Delaunay triangle (this conversion is
#' not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles whose union constitutes the
#' convex hull of \code{Yp} points). Loops are not allowed so arcs are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}) more on the CS-PCDs.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates which serves as a center in the interior of each Delaunay
#' triangle, default for \eqn{M=(1,1,1)} which is the center of mass of each triangle.
#' @param asp A \code{numeric} value, giving the aspect ratio \eqn{y/x} (default is \code{NA}), see the official help page for \code{asp} by typing "\code{? asp}"
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges (default=\code{NULL} for both)
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of the CS-PCD whose vertices are the points in data set \code{Xp} and the Delaunay
#' triangles based on \code{Yp} points
#'
#' @seealso \code{\link{plotCSarcs.tri}}, \code{\link{plotASarcs}}, and \code{\link{plotPEarcs}}
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
#' t<-1.5  #try also t<-2
#'
#' plotCSarcs(Xp,Yp,t,M,xlab="",ylab="")
#' }
#'
#' @export
plotCSarcs <- function(Xp,Yp,t,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
{
  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (nrow(Yp)==3)
  {
    plotCSarcs.tri(Xp,Yp,t,M,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    arcsCS<-arcsCS(Xp,Yp,t,M)
    S<-arcsCS$S
    E<-arcsCS$E

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    Xch<-Xin.convex.hullY(Xp,Yp)

    if (is.null(main))
    {Mvec= paste(M, collapse=",")
    main=paste("Arcs of CS-PCD\n with t = ",t," and M = (",Mvec,")",sep="")
    }

    Xlim<-xlim; Ylim<-ylim
    if (is.null(xlim))
    {xlim<-range(Yp[,1],Xp[,1])
    xr<-xlim[2]-xlim[1]
    xlim<-xlim+xr*c(-.05,.05)
    }
    if (is.null(ylim))
    {ylim<-range(Yp[,2],Xp[,2])
    yr<-ylim[2]-ylim[1]
    ylim<-ylim+yr*c(-.05,.05)
    }

    plot(rbind(Xp),asp=asp,main=main, xlab=xlab, ylab=ylab,xlim=xlim,ylim=ylim,pch=".",cex=3,...)
    interp::plot.triSht(DTmesh, add=TRUE, do.points = TRUE)
    if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}
  }
} #end of the function
#'

#################################################################

#' @title The plot of the Central Similarity (CS) Proximity Regions for a 2D data set - multiple triangle case
#'
#' @description Plots the points in and outside of the Delaunay triangles based on \code{Yp} points which partition
#' the convex hull of \code{Yp} points and also plots the CS proximity regions
#' for \code{Xp} points and the Delaunay triangles based on \code{Yp} points.
#'
#' CS proximity regions are constructed with respect to the Delaunay triangles with the expansion parameter \eqn{t>0}.
#'
#' Edge regions in each triangle is based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)} which is the center of mass of the triangle).
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}) more on the CS proximity regions.
#' Also see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds}) for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points for which CS proximity regions are constructed.
#' @param Yp A set of 2D points which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri} or the circumcenter of \code{tri}.
#' @param asp A \code{numeric} value, giving the aspect ratio \eqn{y/x} (default is \code{NA}), see the official help page for \code{asp} by
#' typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes, respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2, giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the \code{Xp} points, Delaunay triangles based on \code{Yp} and also the CS proximity regions
#' for \code{Xp} points inside the convex hull of \code{Yp} points
#'
#' @seealso \code{\link{plotCSregs.tri}}, \code{\link{plotASregs}} and \code{\link{plotPEregs}}
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
#' tau<-1.5  #try also tau<-2
#'
#' plotCSregs(Xp,Yp,tau,M,xlab="",ylab="")
#' }
#'
#' @export plotCSregs
plotCSregs <- function(Xp,Yp,t,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
{
  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('Xp andYp must be numeric')}

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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (nrow(Yp)==3)
  {
    plotCSregs.tri(Xp,Yp,t,M,asp,main,xlab,ylab,xlim,ylim)
  } else
  { if (!is.point(M,3) || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)  #the Xp points inside the convex hull of Yp points

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)  #number of Delaunay triangles
    nx2<-nrow(Xch)  #number of Xp points inside the convex hull of Yp points

    if (nx2>=1)
    {
      i.tr<-rep(0,nx2)  #the vector of indices for the triangles that contain the Xp points
      for (i1 in 1:nx2)
        for (j1 in 1:nt)
        {
          Tri<-Yp[DTr[j1,],]
          if (in.triangle(Xch[i1,],Tri,boundary=TRUE)$in.tri )
            i.tr[i1]<-j1
        }
    }

    Xlim<-xlim; Ylim<-ylim
    if (is.null(xlim))
    {xlim<-range(Yp[,1],Xp[,1])
    xr<-xlim[2]-xlim[1]
    Xlim<-xlim+xr*c(-.05,.05)
    }
    if (is.null(ylim))
    {ylim<-range(Yp[,2],Xp[,2])
    yr<-ylim[2]-ylim[1]
    Ylim<-ylim+yr*c(-.05,.05)
    }

    if (is.null(main))
    {Mvec= paste(M, collapse=",")
    main=paste("CS Proximity Regions\n with t = ",t," and M = (",Mvec,")",sep="")
    }

    plot(rbind(Xp),main=main, xlab=xlab, ylab=ylab,xlim=Xlim,ylim=Ylim,pch=".",cex=3,...)

    for (i in 1:nt)
    {
      Tri<-Yp[DTr[i,],]  #vertices of the ith triangle
      tri<-as.basic.tri(Tri)$tri #convert the triangle Tri into an nonscaled basic triangle, see as.basic.tri help page

      polygon(tri,lty=2)
      if (nx2>=1)
      {
        Xtri<-matrix(Xch[i.tr==i,],ncol=2)  #Xp points inside triangle i
        ni<-nrow(Xtri)
        if (ni>=1)
        {
          ################
          for (j in 1:ni)
          {
            P1<-Xtri[j,]
            RE<-rel.edge.tri(P1,tri,M)$re

            pr<-NCStri(P1,tri,t,M,re=RE)
            polygon(pr,border="blue")
          }
          ################
        }
      }
    }
  }
} #end of the function
#'

#################################################################
#Domination number functions for NCS
#################################################################

#' @title The indicator for the presence of an arc from a point in set \code{S} to the point \code{p} for Central Similarity
#' Proximity Catch Digraphs (CS-PCDs) - standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p} in \eqn{N_{CS}(x,t)} for some \eqn{x} in \code{S}\eqn{)}, that is, returns 1 if \code{p} is in \eqn{\cup_{x in S} N_{CS}(x,t)},
#' returns 0 otherwise, CS proximity region is constructed with respect to the standard equilateral triangle
#' \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with the expansion parameter \eqn{t>0} and edge regions are based
#' on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e} (which is equivalent to circumcenter of \eqn{T_e}).
#'
#' Edges of \eqn{T_e}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as edges 3, 1, and 2, respectively.
#' If \code{p} is not in \code{S} and either \code{p} or all points in \code{S} are outside \eqn{T_e}, it returns 0,
#' but if \code{p} is in \code{S}, then it always returns 1 regardless of its location (i.e., loops are allowed).
#'
#' See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points. Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param p A 2D point. Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region in the
#' standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)}
#' i.e., the center of mass of \eqn{T_e}.
#'
#' @return \eqn{I(}\code{p} is in \eqn{\cup_{x in S} N_{CS}(x,t))}, that is, returns 1 if \code{p} is in \code{S} or inside \eqn{N_{CS}(x,t)} for at least
#' one \eqn{x} in \code{S}, returns 0 otherwise. CS proximity region is constructed with respect to the standard
#' equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with \code{M}-edge regions.
#'
#' @seealso \code{\link{IarcCSset2pnt.tri}}, \code{\link{IarcCSstd.tri}}, \code{\link{IarcCStri}}, and \code{\link{IarcPEset2pnt.std.tri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-.5
#'
#' S<-rbind(Xp[1,],Xp[2,])  #try also S<-c(.5,.5)
#' IarcCSset2pnt.std.tri(S,Xp[3,],t,M)
#' IarcCSset2pnt.std.tri(S,Xp[3,],t=1,M)
#' IarcCSset2pnt.std.tri(S,Xp[3,],t=1.5,M)
#'
#' S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
#' IarcCSset2pnt.std.tri(S,Xp[3,],t,M)
#' }
#'
#' @export IarcCSset2pnt.std.tri
IarcCSset2pnt.std.tri <- function(S,p,t,M=c(1,1,1))
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.numeric(as.matrix(S)))
  {stop('S must be a matrix of numeric values')}

  if (is.point(S))
  { S<-matrix(S,ncol=2)
  } else
  {S<-as.matrix(S)
  if (ncol(S)!=2 )
  {stop('S must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  k<-nrow(S);
  dom<-0; i<-1;
  while (dom ==0 && i<= k)
  {
    if (IarcCSstd.tri(S[i,],p,t,M)==1)
    {dom<-1};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point in set \code{S} to the point \code{p} for
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) - one triangle case
#'
#' @description Returns I(\code{p} in \eqn{N_{CS}(x,t)} for some \eqn{x} in \code{S}),
#' that is, returns 1 if \code{p} in \eqn{\cup_{x in S} N_{CS}(x,t)},
#' returns 0 otherwise.
#'
#' CS proximity region is constructed with respect to the triangle \code{tri} with
#' the expansion parameter \eqn{t>0} and edge regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' Edges of \code{tri}\eqn{=T(A,B,C)}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as edges 3, 1, and 2, respectively.
#' If \code{p} is not in \code{S} and either \code{p} or all points in \code{S} are outside \code{tri}, it returns 0,
#' but if \code{p} is in \code{S}, then it always returns 1 regardless of its location (i.e., loops are allowed).
#'
#' @param S A set of 2D points. Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param p A 2D point. Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region
#' constructed in the triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' @return I(\code{p} is in \eqn{\cup_{x in S} N_{CS}(x,t)}), that is, returns 1 if \code{p} is in \code{S} or inside \eqn{N_{CS}(x,t)} for at least
#' one \eqn{x} in \code{S}, returns 0 otherwise where CS proximity region is constructed with respect to the triangle \code{tri}
#'
#' @seealso \code{\link{IarcCSset2pnt.std.tri}}, \code{\link{IarcCStri}}, \code{\link{IarcCSstd.tri}},
#' \code{\link{IarcASset2pnt.tri}}, and \code{\link{IarcPEset2pnt.tri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' tau<-.5
#'
#' IarcCSset2pnt.tri(S,Xp[3,],Tr,tau,M)
#' IarcCSset2pnt.tri(S,Xp[3,],Tr,t=1,M)
#' IarcCSset2pnt.tri(S,Xp[3,],Tr,t=1.5,M)
#'
#' S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
#' IarcCSset2pnt.tri(S,Xp[3,],Tr,tau,M)
#' }
#'
#' @export IarcCSset2pnt.tri
IarcCSset2pnt.tri <- function(S,p,tri,t,M=c(1,1,1))
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.numeric(as.matrix(S)))
  {stop('S must be a matrix of numeric values')}

  if (is.point(S))
  { S<-matrix(S,ncol=2)
  } else
  {S<-as.matrix(S)
  if (ncol(S)!=2 )
  {stop('S must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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

  k<-nrow(S);
  dom<-0; i<-1;
  while (dom ==0 && i<= k)
  {
    if (IarcCStri(S[i,],p,tri,t,M)==1)
    {dom<-1};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the set of points \code{S} being a dominating set or not for Central Similarity Proximity
#' Catch Digraphs (CS-PCDs) - standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{S} a dominating set of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the data set \code{Xp}), that is,
#' returns 1 if \code{S} is a dominating set of CS-PCD, returns 0 otherwise.
#'
#' CS proximity region is constructed
#' with respect to the standard equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with
#' expansion parameter \eqn{t>0} and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e} (which is equivalent to the circumcenter of \eqn{T_e}).
#'
#' Edges of \eqn{T_e}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as 3, 1, and 2, respectively.
#'
#' See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points which is to be tested for being a dominating set for the CS-PCDs.
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region in the
#' standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return \eqn{I(}\code{S} a dominating set of the CS-PCD\eqn{)}, that is, returns 1 if \code{S} is a dominating set of CS-PCD,
#' returns 0 otherwise, where CS proximity region is constructed in the standard equilateral triangle \eqn{T_e}
#'
#' @seealso \code{\link{Idom.setCStri}} and \code{\link{Idom.setPEstd.tri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-.5
#'
#' S<-rbind(Xp[1,],Xp[2,])
#' Idom.setCSstd.tri(S,Xp,t,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' Idom.setCSstd.tri(S,Xp,t,M)
#' }
#'
#' @export Idom.setCSstd.tri
Idom.setCSstd.tri <- function(S,Xp,t,M=c(1,1,1))
{
  if (!is.numeric(as.matrix(S)) || !is.numeric(as.matrix(Xp)))
  {stop('Both arguments must be numeric')}

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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  k<-nrow(S);
  n<-nrow(Xp);

  dom<-1; i<-1;
  while (dom ==1 && i<= n)
  {
    if (IarcCSset2pnt.std.tri(S,Xp[i,],t,M)==0)  #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the set of points \code{S} being a dominating set or not for Central Similarity Proximity
#' Catch Digraphs (CS-PCDs) - one triangle case
#'
#' @description Returns \eqn{I(}\code{S} a dominating set of CS-PCD whose vertices are the data set \code{Xp}\eqn{)}, that is,
#' returns 1 if \code{S} is a dominating set of CS-PCD, returns 0 otherwise.
#'
#' CS proximity region is constructed with
#' respect to the triangle \code{tri} with the expansion parameter \eqn{t>0} and edge regions are based
#' on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}; default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' The triangle \code{tri}\eqn{=T(A,B,C)} has edges \eqn{AB}, \eqn{BC}, \eqn{AC} which are also labeled as edges 3, 1, and 2, respectively.
#'
#' See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points which is to be tested for being a dominating set for the CS-PCDs.
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region
#' constructed in the triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' @return \eqn{I(}\code{S} a dominating set of the CS-PCD\eqn{)}, that is, returns 1 if \code{S} is a dominating set of CS-PCD whose
#' vertices are the data points in \code{Xp}; returns 0 otherwise, where CS proximity region is constructed in
#' the triangle \code{tri}
#'
#' @seealso \code{\link{Idom.setCSstd.tri}}, \code{\link{Idom.setPEtri}} and \code{\link{Idom.setAStri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' tau<-.5
#' S<-rbind(Xp[1,],Xp[2,])
#' Idom.setCStri(S,Xp,Tr,tau,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' Idom.setCStri(S,Xp,Tr,tau,M)
#' }
#'
#' @export Idom.setCStri
Idom.setCStri <- function(S,Xp,tri,t,M=c(1,1,1))
{
  if (!is.numeric(as.matrix(S)) || !is.numeric(as.matrix(Xp)))
  {stop('Both arguments must be numeric')}

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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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

  if (isTRUE(all.equal(M,circumcenter.tri(tri)))==FALSE & in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  k<-nrow(S);
  n<-nrow(Xp);

  dom<-1; i<-1;
  while (dom ==1 && i<= n)
  {
    if (IarcCSset2pnt.tri(S,Xp[i,],tri,t,M)==0)  #this is where tri is used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for \code{k} being an upper bound for the domination number of Central Similarity Proximity
#' Catch Digraph (CS-PCD) by the exact algorithm - standard equilateral triangle case
#'
#' @description Returns \eqn{I(}domination number of CS-PCD is less than or equal to \code{k}\eqn{)} where the vertices of the CS-PCD are the data points \code{Xp},
#' that is, returns 1 if the domination number of CS-PCD is less than the prespecified value \code{k}, returns 0
#' otherwise. It also provides the vertices (i.e., data points) in a dominating set of size \code{k} of CS-PCD.
#'
#' CS proximity region is constructed with respect to the standard equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with
#' expansion parameter \eqn{t>0} and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e} (which is equivalent to the circumcenter of \eqn{T_e}).
#'
#' Edges of \eqn{T_e}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as 3, 1, and 2, respectively.
#' Loops are allowed in the digraph.
#' It takes a long time for large number of vertices (i.e., large number of row numbers).
#'
#' See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points which constitute the vertices of CS-PCD.
#' @param k A positive integer representing an upper bound for the domination number of CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region in the
#' standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard equilateral triangle \eqn{T_e}; default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return A \code{list} with two elements
#' \item{domUB}{The upper bound \code{k} (to be checked) for the domination number of CS-PCD. It is prespecified
#' as \code{k} in the function arguments.}
#' \item{Idom.num.up.bnd}{The indicator for the upper bound for domination number of CS-PCD being the
#' specified value \code{k} or not. It returns 1 if the upper bound is \code{k}, and 0 otherwise.}
#' \item{ind.domset}{The vertices (i.e., data points) in the dominating set of size \code{k} if it exists,
#' otherwise it is \code{NULL}.}
#'
#' @seealso \code{\link{Idom.numCSup.bnd.tri}}, \code{\link{Idom.num.up.bnd}}, \code{\link{Idom.numASup.bnd.tri}},
#' and \code{\link{dom.num.exact}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' t<-.5
#'
#' Idom.numCSup.bnd.std.tri(Xp,1,t,M)
#'
#' for (k in 1:n)
#'   print(c(k,Idom.numCSup.bnd.std.tri(Xp,k,t,M)$Idom.num.up.bnd))
#'   print(c(k,Idom.numCSup.bnd.std.tri(Xp,k,t,M)$domUB))
#' }
#'
#' @export Idom.numCSup.bnd.std.tri
Idom.numCSup.bnd.std.tri <- function(Xp,k,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp);
  xc<-combinat::combn(1:n,k); N1<-choose(n,k);
  xc<-matrix(xc,ncol=N1)
  dom<-0; j<-1; domset<-c();
  while (j<=N1 && dom==0)
  {
    if (Idom.setCSstd.tri(Xp[xc[,j],],Xp,t,M)==1)  #this is where std equilateral triangle Te is implicitly used
    {dom<-1; domset<-Xp[xc[,j],];}
    j<-j+1;
  }

  list(domUB=k, #upper bound for the domination number of CS-PCD
       Idom.num.up.bnd=dom, #indicator that domination number <=k
       domset=domset #a dominating set of size k (if exists)
  )
} #end of the function
#'

#################################################################

#' @title Indicator for an upper bound for the domination number of Central Similarity Proximity Catch Digraph
#' (CS-PCD) by the exact algorithm - one triangle case
#'
#' @description Returns \eqn{I(}domination number of CS-PCD is less than or equal to \code{k}\eqn{)} where the vertices of the CS-PCD are the data points \code{Xp},
#' that is, returns 1 if the domination number of CS-PCD is less than the prespecified value \code{k}, returns 0
#' otherwise. It also provides the vertices (i.e., data points) in a dominating set of size \code{k} of CS-PCD.
#'
#' CS proximity region is constructed with respect to the triangle \code{tri}\eqn{=T(A,B,C)} with
#' expansion parameter \eqn{t>0} and edge regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#'
#' Edges of \code{tri}, \eqn{AB}, \eqn{BC}, \eqn{AC}, are also labeled as 3, 1, and 2, respectively.
#' Loops are allowed in the digraph.
#'
#' See also (\insertCite{ceyhan:mcap2012;textual}{pcds}).
#'
#' Caveat: It takes a long time for large number of vertices (i.e., large number of row numbers).
#'
#' @param Xp A set of 2D points which constitute the vertices of CS-PCD.
#' @param k A positive integer to be tested for an upper bound for the domination number of CS-PCDs.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region in the
#' triangle \code{tri}.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates which serves as a
#' center in the interior of the triangle \code{tri}; default is \eqn{M=(1,1,1)}, i.e.
#' the center of mass of \code{tri}.
#'
#' @return A \code{list} with two elements
#' \item{domUB}{The upper bound \code{k} (to be checked) for the domination number of CS-PCD. It is prespecified
#' as \code{k} in the function arguments.}
#' \item{Idom.num.up.bnd}{The indicator for the upper bound for domination number of CS-PCD being the
#' specified value \code{k} or not. It returns 1 if the upper bound is \code{k}, and 0 otherwise.}
#' \item{ind.domset}{The vertices (i.e., data points) in the dominating set of size \code{k} if it exists,
#' otherwise it is \code{NULL}.}
#'
#' @seealso \code{\link{Idom.numCSup.bnd.std.tri}}, \code{\link{Idom.num.up.bnd}}, \code{\link{Idom.numASup.bnd.tri}},
#' and \code{\link{dom.num.exact}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' t<-.5
#'
#' Idom.numCSup.bnd.tri(Xp,1,Tr,t,M)
#'
#' for (k in 1:n)
#'   print(c(k,Idom.numCSup.bnd.tri(Xp,k,Tr,t,M)))
#' }
#'
#' @export Idom.numCSup.bnd.tri
Idom.numCSup.bnd.tri <- function(Xp,k,tri,t,M=c(1,1,1))
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

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

  n<-nrow(Xp);
  xc<-combinat::combn(1:n,k); N1<-choose(n,k);
  xc<-matrix(xc,ncol=N1)
  dom<-0; j<-1; domset<-c();
  while (j<=N1 && dom==0)
  {
    if (Idom.setCStri(Xp[xc[,j],],Xp,tri,t,M)==1)  #this is where triangle tri is used
    {dom<-1; domset<-Xp[xc[,j],];}
    j<-j+1;
  }

  list(domUB=k, #upper bound for the domination number of CS-PCD
       Idom.num.up.bnd=dom, #indicator that domination number <=k
       domset=domset #a dominating set of size k (if exists)
  )
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for Central Similarity Proximity Catch Digraphs
#' (CS-PCDs) - standard equilateral triangle case with \eqn{t=1}
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp} in the standard equilateral
#' triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}, that is, returns 1 if \code{p} is a dominating point of CS-PCD,
#' returns 0 otherwise.
#'
#' Point, \code{p}, is in the edge region of edge re (default is \code{NULL}); vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise in \eqn{T_e}, and the opposite edges are labeled with label of the vertices
#' (that is, edge numbering is \code{1,2}, and \code{3} for edges \eqn{AB}, \eqn{BC}, and \eqn{AC}).
#'
#' CS proximity region is constructed with respect to \eqn{T_e} with expansion parameter \eqn{t=1}
#' and edge regions are based on center of mass \eqn{CM=(1/2,\sqrt{3}/6)}.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}), so by default this
#' function checks whether the point \code{p} would be a dominating point if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point or not of the CS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the CS-PCD.
#' @param re The index of the edge region in \eqn{T_e} containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise.
#'
#' @seealso \code{\link{Idom.num1CSstd.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' Te<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num1CSt1std.tri(Xp[3,],Xp)
#'
#' Idom.num1CSt1std.tri(c(1,2),c(1,2))
#' Idom.num1CSt1std.tri(c(1,2),c(1,2),ch.data.pnt = TRUE)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1CSt1std.tri(Xp[i,],Xp))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp)
#' L<-Te; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE);
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(rbind(Xp[ind.gam1,]),pch=4,col=2)
#' #rbind is to insert the points correctly if there is only one dominating point
#'
#' txt<-rbind(Te,CM)
#' xc<-txt[,1]+c(-.02,.02,.01,.05)
#' yc<-txt[,2]+c(.02,.02,.03,.02)
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
Idom.num1CSt1std.tri <- function(p,Xp,re=NULL,ch.data.pnt=FALSE)
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

  if (length(Xp)==2 && isTRUE(all.equal(matrix(p,ncol=2),matrix(Xp,ncol=2))) )
  {dom<-1; return(dom); stop}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C)
  if (!in.triangle(p,Te, boundary = TRUE)$in.tri)
  {dom<-0; return(dom); stop}

  if (is.null(re))
  {re<-rel.edge.std.triCM(p)$re #edge region for p
  } else
  {  if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')}}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('p is not a data point in Xp')}
  }

  dom<-0
  c.e<-cl2edges.std.tri(Xp)$ext
  c.e<-c.e[-re,]
  if (IarcCSt1.std.tri(p,c.e[1,])==1 && IarcCSt1.std.tri(p,c.e[2,])==1)  #this is where std equilateral triangle Te is implicitly used
    dom<-1;
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for Central Similarity Proximity Catch Digraphs
#' (CS-PCDs) - standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp} in the standard equilateral
#' triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}, that is, returns 1 if \code{p} is a dominating point of CS-PCD,
#' returns 0 otherwise.
#'
#' CS proximity region is constructed with respect to \eqn{T_e} with expansion parameter \eqn{t>0}
#' and edge regions are based on center of mass \eqn{CM=(1/2,\sqrt{3}/6)}.
#'
#' \code{ch.data.pnt} is for checking whether point \code{p} is a data point in \code{Xp} or not (default is \code{FALSE}), so by default this
#' function checks whether the point \code{p} would be a dominating point if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point or not of the CS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param ch.data.pnt A logical argument for checking whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the CS-PCD\eqn{)} where the vertices of the CS-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point, returns 0 otherwise
#'
#' @seealso \code{\link{Idom.num1CSt1std.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' CM<-(A+B+C)/3
#' Te<-rbind(A,B,C);
#' t<-1.5
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num1CSstd.tri(Xp[3,],Xp,t)
#' Idom.num1CSstd.tri(c(1,2),c(1,2),t)
#' Idom.num1CSstd.tri(c(1,2),c(1,2),t,ch.data.pnt = TRUE)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1CSstd.tri(Xp[i,],Xp,t))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Te,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp)
#' L<-Te; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE);
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(rbind(Xp[ind.gam1,]),pch=4,col=2)
#' #rbind is to insert the points correctly if there is only one dominating point
#'
#' txt<-rbind(Te,CM)
#' xc<-txt[,1]+c(-.02,.02,.01,.05)
#' yc<-txt[,2]+c(.02,.02,.03,.02)
#' txt.str<-c("A","B","C","CM")
#' text(xc,yc,txt.str)
#'
#' Idom.num1CSstd.tri(c(1,2),Xp,t,ch.data.pnt = FALSE)
#' #gives an error if ch.data.pnt = TRUE message since p is not a data point
#' }
#'
#' @export Idom.num1CSstd.tri
Idom.num1CSstd.tri <- function(p,Xp,t,ch.data.pnt=FALSE)
{
  if (!is.point(p))
  {stop('p must be a numeric point of dimension 2')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('second argument must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('point, p, is not a data point in Xp')}
  }

  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n && dom==1)
  {if (IarcCSstd.tri(p,Xp[i,],t)==0)  #this is where std equilateral triangle Te is implicitly used
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

# funsCSGamTe
#'
#' @title The function \code{gammakCSstd.tri} is for \eqn{k} (\eqn{k=2,3,4,5}) points constituting a dominating set for Central Similarity
#' Proximity Catch Digraphs (CS-PCDs) - standard equilateral triangle case
#'
#' @description
#' Four indicator functions: \code{Idom.num2CSstd.tri}, \code{Idom.num3CSstd.tri}, \code{Idom.num4CSstd.tri}, \code{Idom.num5CSstd.tri} and \code{Idom.num6CSstd.tri}.
#'
#' The function \code{gammakCSstd.tri} returns I(\{\code{p1},...,\code{pk}\} is a dominating set of the CS-PCD)
#' where vertices of CS-PCD are the 2D data set \code{Xp}, that is, returns 1 if \{\code{p1},...,\code{pk}\}
#' is a dominating set of CS-PCD, returns 0 otherwise for \eqn{k=2,3,4,5,6}.
#'
#' CS proximity region is constructed with respect to \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' with expansion parameter \eqn{t>0} and edge regions are based on center of mass \eqn{CM=(1/2,\sqrt{3}/6)}.
#'
#' \code{ch.data.pnts} is for checking whether points \code{p1},...,\code{pk} are data points in \code{Xp} or not
#' (default is \code{FALSE}), so by default this function checks whether the points \code{p1},...,\code{pk} would be a
#' dominating set if they actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p1,p2,p3,p4,p5,p6 The points \{\eqn{p1,\ldots,pk}\} are \eqn{k} 2D points
#' (for \eqn{k=2,3,4,5,6}) to be tested for constituting a dominating set of the CS-PCD.
#' @param Xp A set of 2D points which constitutes the vertices of the CS-PCD.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param ch.data.pnts A logical argument for checking whether points \{\eqn{p1,\ldots,pk}\} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return The function \code{gammakCSstd.tri} returns \{\code{p1},...,\code{pk}\} is a dominating set of the CS-PCD) where
#' vertices of the CS-PCD are the 2D data set \code{Xp}), that is, returns 1 if \{\code{p1},...,\code{pk}\}
#' is a dominating set of CS-PCD, returns 0 otherwise.
#'
#' @name funsCSGamTe
NULL
#'
#' @seealso \code{\link{Idom.num1CSstd.tri}}, \code{\link{Idom.num2PEtri}} and \code{\link{Idom.num2PEtetra}}
#'
#' @rdname funsCSGamTe
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' #Examples for Idom.num2CSstd.tri
#' t<-1.5
#' n<-10 #try also 10, 20 (it may take longer for larger n)
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num2CSstd.tri(Xp[1,],Xp[2,],Xp,t)
#' Idom.num2CSstd.tri(c(.2,.2),Xp[2,],Xp,t)
#'
#' ind.gam2<-vector()
#' for (i in 1:(n-1))
#'  for (j in (i+1):n)
#'  {if (Idom.num2CSstd.tri(Xp[i,],Xp[j,],Xp,t)==1)
#'   ind.gam2<-rbind(ind.gam2,c(i,j))}
#'
#' ind.gam2
#' }
#'
#' @export
Idom.num2CSstd.tri <- function(p1,p2,Xp,t,ch.data.pnts=FALSE)
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

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points are data points in Xp')}
  }

  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i<=n && dom==1)
  {
    if (IarcCSstd.tri(p1,Xp[i,],t)==0 && IarcCSstd.tri(p2,Xp[i,],t)==0)  #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'
#' @rdname funsCSGamTe
#'
#' @examples
#' \dontrun{
#' #Examples for Idom.num3CSstd.tri
#' t<-1.5
#' n<-10 #try also 10, 20 (it may take longer for larger n)
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num3CSstd.tri(Xp[1,],Xp[2,],Xp[3,],Xp,t)
#'
#' ind.gam3<-vector()
#' for (i in 1:(n-2))
#'  for (j in (i+1):(n-1))
#'    for (k in (j+1):n)
#'    {if (Idom.num3CSstd.tri(Xp[i,],Xp[j,],Xp[k,],Xp,t)==1)
#'     ind.gam3<-rbind(ind.gam3,c(i,j,k))}
#'
#' ind.gam3
#' }
#'
#' @export
Idom.num3CSstd.tri <- function(p1,p2,p3,Xp,t,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) || !is.point(p3) )
  {stop('p1, p2, and p3 must all be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp) || !is.in.data(p3,Xp))
    {stop('not all points are data points in Xp')}
  }

  n<-nrow(Xp);
  dom<-1; i<-1;
  while (i<=n && dom==1)
  {
    if (IarcCSstd.tri(p1,Xp[i,],t)==0 && IarcCSstd.tri(p2,Xp[i,],t)==0 && IarcCSstd.tri(p3,Xp[i,],t)==0)  #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'
#' @rdname funsCSGamTe
#'
#' @examples
#' \dontrun{
#' #Examples for Idom.num4CSstd.tri
#' t<-1.5
#' n<-10 #try also 10, 20 (it may take longer for larger n)
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num4CSstd.tri(Xp[1,],Xp[2,],Xp[3,],Xp[4,],Xp,t)
#'
#' ind.gam4<-vector()
#' for (i in 1:(n-3))
#'  for (j in (i+1):(n-2))
#'    for (k in (j+1):(n-1))
#'      for (l in (k+1):n)
#'      {if (Idom.num4CSstd.tri(Xp[i,],Xp[j,],Xp[k,],Xp[l,],Xp,t)==1)
#'       ind.gam4<-rbind(ind.gam4,c(i,j,k,l))}
#'
#' ind.gam4
#'
#' Idom.num4CSstd.tri(c(.2,.2),Xp[2,],Xp[3,],Xp[4,],Xp,t,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE since not all points are data points in Xp
#' }
#'
#' @export
Idom.num4CSstd.tri <- function(p1,p2,p3,p4,Xp,t,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) || !is.point(p3) || !is.point(p4) )
  {stop('p1, p2, p3 and p4 must all be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp) || !is.in.data(p3,Xp) || !is.in.data(p4,Xp))
    {stop('not all points are data points in Xp')}
  }

  n<-nrow(Xp);
  dom<-1; i<-1;
  while (i<=n && dom==1)
  {
    if (IarcCSstd.tri(p1,Xp[i,],t)==0 && IarcCSstd.tri(p2,Xp[i,],t)==0 && IarcCSstd.tri(p3,Xp[i,],t)==0
        && IarcCSstd.tri(p4,Xp[i,],t)==0)
      #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'
#' @rdname funsCSGamTe
#'
#' @examples
#' \dontrun{
#' #Examples for Idom.num5CSstd.tri
#' t<-1.5
#' n<-10 #try also 10, 20 (it may take longer for larger n)
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num5CSstd.tri(Xp[1,],Xp[2,],Xp[3,],Xp[4,],Xp[5,],Xp,t)
#'
#' ind.gam5<-vector()
#' for (i1 in 1:(n-4))
#'  for (i2 in (i1+1):(n-3))
#'    for (i3 in (i2+1):(n-2))
#'      for (i4 in (i3+1):(n-1))
#'        for (i5 in (i4+1):n)
#'        {if (Idom.num5CSstd.tri(Xp[i1,],Xp[i2,],Xp[i3,],Xp[i4,],Xp[i5,],Xp,t)==1)
#'         ind.gam5<-rbind(ind.gam5,c(i1,i2,i3,i4,i5))}
#'
#' ind.gam5
#'
#' Idom.num5CSstd.tri(c(.2,.2),Xp[2,],Xp[3,],Xp[4,],Xp[5,],Xp,t,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE since not all points are data points in Xp
#' }
#'
#' @export
Idom.num5CSstd.tri <- function(p1,p2,p3,p4,p5,Xp,t,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) || !is.point(p3) || !is.point(p4) || !is.point(p5) )
  {stop('p1, p2, p3, p4 and p5 must all be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp) || !is.in.data(p3,Xp) || !is.in.data(p4,Xp) || !is.in.data(p5,Xp))
    {stop('not all points are data points in Xp')}
  }

  n<-nrow(Xp);
  dom<-1; i<-1;
  while (i<=n && dom==1)
  {
    if (IarcCSstd.tri(p1,Xp[i,],t)==0 && IarcCSstd.tri(p2,Xp[i,],t)==0 && IarcCSstd.tri(p3,Xp[i,],t)==0 &&
        IarcCSstd.tri(p4,Xp[i,],t)==0 && IarcCSstd.tri(p5,Xp[i,],t)==0)  #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'
#' @rdname funsCSGamTe
#'
#' @examples
#' \dontrun{
#' #Examples for Idom.num6CSstd.tri
#' t<-1.5
#' n<-10 #try also 10, 20 (it may take longer for larger n)
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Idom.num6CSstd.tri(Xp[1,],Xp[2,],Xp[3,],Xp[4,],Xp[5,],Xp[6,],Xp,t)
#'
#' ind.gam6<-vector()
#' for (i1 in 1:(n-5))
#'  for (i2 in (i1+1):(n-4))
#'    for (i3 in (i2+1):(n-3))
#'      for (i4 in (i3+1):(n-2))
#'        for (i5 in (i4+1):(n-1))
#'          for (i6 in (i5+1):n)
#'          {if (Idom.num6CSstd.tri(Xp[i1,],Xp[i2,],Xp[i3,],Xp[i4,],Xp[i5,],Xp[i6,],Xp,t)==1)
#'           ind.gam6<-rbind(ind.gam6,c(i1,i2,i3,i4,i5,i6))}
#'
#' ind.gam6
#'
#' Idom.num6CSstd.tri(c(.2,.2),Xp[2,],Xp[3,],Xp[4,],Xp[5,],Xp[6,],Xp,t,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE since not all points are data points in Xp
#' }
#'
#' @export
Idom.num6CSstd.tri <- function(p1,p2,p3,p4,p5,p6,Xp,t,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2) || !is.point(p3) || !is.point(p4) || !is.point(p5) || !is.point(p6) )
  {stop('p1, p2, p3, p4, p5 and p6 must all be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(t,1) || t<=0)
  {stop('r must be a scalar greater than 0')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp) || !is.in.data(p3,Xp) || !is.in.data(p4,Xp)
        || !is.in.data(p5,Xp) || !is.in.data(p6,Xp))
    {stop('not all points are data points in Xp')}
  }

  n<-nrow(Xp);
  dom<-1; i<-1;
  while (i<=n && dom==1)
  {
    if (IarcCSstd.tri(p1,Xp[i,],t)==0 && IarcCSstd.tri(p2,Xp[i,],t)==0 && IarcCSstd.tri(p3,Xp[i,],t)==0 &&
        IarcCSstd.tri(p4,Xp[i,],t)==0 && IarcCSstd.tri(p5,Xp[i,],t)==0 && IarcCSstd.tri(p6,Xp[i,],t)==0)  #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title  An alternative to the function \code{\link{IarcCStri}} which yields the indicator
#' for the presence of an arc from one point to another
#' for Central Similarity Proximity Catch Digraphs (CS-PCDs)
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for points \code{p1} and \code{p2}, that is, returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)},
#' returns 0 otherwise, where \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x} with the expansion parameter \eqn{t>0}.
#'
#' CS proximity region is constructed with respect to the triangle \code{tri} and edge regions are based on the
#' center of mass, \eqn{CM}. \code{re} is the index of the \eqn{CM}-edge region \code{p} resides, with default=\code{NULL} but must be provided as
#' vertices \eqn{(y_1,y_2,y_3)} for \eqn{re=3} as rbind(y2,y3,y1) for \eqn{re=1} and as rbind(y1,y3,y2) for \eqn{re=2} for triangle \eqn{T(y_1,y_2,y_3)}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are outside \code{tri}, it returns 0,
#' but if they are identical, then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}).
#'
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines whether \code{p2} is inside the CS proximity region of
#' \code{p1} or not.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param re Index of the \eqn{CM}-edge region containing the point \code{p},
#' either \code{1,2,3} or \code{NULL}, default=\code{NULL} but
#' must be provided (row-wise) as vertices \eqn{(y_1,y_2,y_3)} for \code{re=3} as \eqn{(y_2,y_3,y_1)} for
#' \code{re=1} and as \eqn{(y_1,y_3,y_2)} for \code{re=2} for triangle \eqn{T(y_1,y_2,y_3)}.
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{CS}(p1,t))} for \code{p1}, that is,
#' returns 1 if \code{p2} is in \eqn{N_{CS}(p1,t)}, returns 0 otherwise.
#'
#' @seealso \code{\link{IarcAStri}}, \code{\link{IarcPEtri}}, \code{\link{IarcCStri}}, and \code{\link{IarcCSstd.tri}}
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
#' t<-1.5
#'
#' P1<-c(.4,.2)
#' P2<-c(1.8,.5)
#' IarcCStri(P1,P2,Tr,t,M=c(1,1,1))
#' IarcCStri.alt(P1,P2,Tr,t)
#'
#' IarcCStri(P2,P1,Tr,t,M=c(1,1,1))
#' IarcCStri.alt(P2,P1,Tr,t)
#'
#' #or try
#' re<-rel.edges.triCM(P1,Tr)$re
#' IarcCStri(P1,P2,Tr,t,M=c(1,1,1),re)
#' IarcCStri.alt(P1,P2,Tr,t,re)
#' }
#'
#' @export
IarcCStri.alt <- function(p1,p2,tri,t,re=NULL)
{
  if (!is.point(p1) || !is.point(p2) )
  {stop('p1 and p2 must both be numeric 2D points')}

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('the triangle is degenerate')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,tri,boundary=TRUE)$in.tri || !in.triangle(p2,tri,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(re))
  {re<-edge.reg.triCM(p1,tri)  #related vertices for edge region for p1
  } else
  {  if (!is.numeric(re) || sum(re==c(1,2,3))!=1)
  {stop('edge index, re, must be 1, 2 or 3')}}

  Arc<-0
  CM<-apply(tri,2,mean);
  Dx <-(p1[2]*CM[1]*re[2,1]-p1[2]*CM[1]*re[1,1]-p1[2]*re[1,1]*re[2,1]+p1[2]*re[1,1]^2-p1[1]*CM[2]*re[2,1]+
          p1[1]*CM[2]*re[1,1]+re[1,2]*p1[1]*re[2,1]-re[1,2]*p1[1]*re[1,1]+re[2,2]*CM[1]*re[1,1]-re[2,2]*re[1,1]^2-
          re[2,1]*re[1,2]*CM[1]+re[2,1]*re[1,2]*re[1,1])/(-CM[2]*re[2,1]+CM[2]*re[1,1]+re[2,1]*re[1,2]+re[2,2]*CM[1]-
                                                            re[2,2]*re[1,1]-re[1,2]*CM[1]);
  Dy<-Line(re[1,],re[2,],Dx)$y;
  Ex<-(-p1[2]*re[2,1]^2+p1[2]*re[1,1]*re[2,1]+p1[2]*CM[1]*re[2,1]-p1[2]*CM[1]*re[1,1]-p1[1]*CM[2]*re[2,1]+
         p1[1]*CM[2]*re[1,1]+re[2,2]*p1[1]*re[2,1]-re[2,2]*p1[1]*re[1,1]-re[2,2]*re[1,1]*re[2,1]+re[2,2]*CM[1]*re[1,1]-
         re[2,1]*re[1,2]*CM[1]+re[2,1]^2*re[1,2])/(-CM[2]*re[2,1]+CM[2]*re[1,1]+re[2,1]*re[1,2]+re[2,2]*CM[1]-
                                                     re[2,2]*re[1,1]-re[1,2]*CM[1]);
  Ey<-Line(re[1,],re[2,],Ex)$y;
  D1<-p1+t*(c(Dx,Dy)-p1);
  E1<-p1+t*(c(Ex,Ey)-p1);

  if (sign(p1[2]-paraline(D1,re[1,],re[3,],p1[1])$y)==sign(p2[2]-paraline(D1,re[1,],re[3,],p2[1])$y) &&
      sign(p1[2]-paraline(E1,re[2,],re[3,],p1[1])$y)==sign(p2[2]-paraline(E1,re[2,],re[3,],p2[1])$y) &&
      sign(p1[2]-Line(D1,E1,p1[1])$y)==sign(p2[2]-Line(D1,E1,p2[1])$y))
    Arc<-1;
  Arc
} #end of the function
#'

#################################################################

#' @title Arc density of Central Similarity Proximity Catch Digraphs (CS-PCDs) - one triangle case
#'
#' @description Returns the arc density of CS-PCD whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the triangle \code{tri}.
#'
#' CS proximity regions is defined with respect to \code{tri} with
#' expansion parameter \eqn{t>0} and edge regions are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}; default is
#' \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' The function also provides arc density standardized by the mean and asymptotic variance of the arc density
#' of CS-PCD for uniform data in the triangle \code{tri}
#' only when \code{M} is the center of mass. For the number of arcs, loops are not allowed.
#'
#' \code{tri.cor} is a logical argument for triangle correction (default is \code{TRUE}), if \code{TRUE}, only the points inside the
#' triangle are considered (i.e., digraph induced by these vertices are considered) in computing the arc density,
#' otherwise all points are considered (for the number of vertices in the denominator of arc density).
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS,ceyhan:test2014;textual}{pcds}) for more on CS-PCDs.
#'
#' @param Xp A set of 2D points which constitute the vertices of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param t A positive real number which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \code{tri}.
#' @param tri.cor A logical argument for computing the arc density for only the points inside the triangle,
#' \code{tri} (default is \code{tri.cor=FALSE}), i.e., if \code{tri.cor=TRUE} only the induced digraph with the vertices
#' inside \code{tri} are considered in the computation of arc density.
#'
#' @return A \code{list} with the elements
#' \item{arc.dens}{Arc density of CS-PCD whose vertices are the 2D numerical data set, \code{Xp};
#' CS proximity regions are defined with respect to the triangle \code{tri} and \code{M}-edge regions}
#' \item{std.arc.dens}{Arc density standardized by the mean and asymptotic variance of the arc
#' density of CS-PCD for uniform data in the triangle \code{tri}.This will only be returned if \code{M} is the center of mass.}
#'
#' @seealso \code{\link{ASarc.dens.tri}}, \code{\link{PEarc.dens.tri}}, and \code{\link{num.arcsCStri}}
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
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' CSarc.dens.tri(Xp,Tr,t=.5,M)
#' CSarc.dens.tri(Xp,Tr,t=.5,M,tri.cor = FALSE)
#' #try also t=1 and t=1.5 above
#' }
#'
#' @export CSarc.dens.tri
CSarc.dens.tri <- function(Xp,tri,t,M=c(1,1,1),tri.cor=FALSE)
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

  if (!is.point(t,1) || t<=0)
  {stop('t must be a scalar greater than 0')}

  if (!is.point(M) && !is.point(M,3) )
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (isTRUE(all.equal(M,circumcenter.tri(tri)))==FALSE & in.triangle(M,tri,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  nx<-nrow(Xp)

  narcs<-num.arcsCStri(Xp,tri,t,M)$num.arcs
  mean.rho<-muCS2D(t)
  var.rho<-asyvarCS2D(t)

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
    {stop('There are not enough Xp points in the triangle, tri, to compute the arc density!')}
    n<-NinTri
  } else
  {
    n<-nx
  }

  rho<-narcs/(n*(n-1))
  res=list(arc.dens=rho)

  CM=apply(tri,2,mean)
  if (isTRUE(all.equal(M,CM))){
    std.rho<-sqrt(n)*(rho-mean.rho)/sqrt(var.rho)
    res=list(
      arc.dens=rho, #arc density
      std.arc.dens=std.rho
    )}

  res
} #end of the function
#'
