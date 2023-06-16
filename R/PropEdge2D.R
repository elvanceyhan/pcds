#PropEdge2D.R;
#Functions for NPE in R^2
#################################################################

#' @title The vertices of the Proportional Edge (PE) Proximity Region
#' in a standard basic triangle
#'
#' @description Returns the vertices of the PE proximity region
#' (which is itself a triangle) for a point in the
#' standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))=}\code{(rv=1,rv=2,rv=3)}.
#'
#' PE proximity region is defined with respect
#' to the standard basic triangle \eqn{T_b}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions based on center \eqn{M=(m_1,m_2)} in
#' Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of the basic
#' triangle \eqn{T_b} or based on the circumcenter of \eqn{T_b};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#'
#' Vertex regions are labeled as \eqn{1,2,3} rowwise for the vertices
#' of the triangle \eqn{T_b}. \code{rv} is the index of the vertex region
#' \code{p} resides, with default=\code{NULL}.
#' If \code{p} is outside of \code{tri},
#' it returns \code{NULL} for the proximity region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point whose PE proximity region is to be computed.
#' @param r A positive real number which serves
#' as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c1,c2 Positive real numbers
#' representing the top vertex in standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0}
#' and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point
#' in barycentric coordinates
#' which serves as a center
#' in the interior of the standard basic triangle \eqn{T_b}
#' or the circumcenter of \eqn{T_b}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#' @param rv Index of the \code{M}-vertex region
#' containing the point \code{p}, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return Vertices of the triangular region
#' which constitutes the PE proximity region with expansion parameter
#' r and center \code{M} for a point \code{p}
#'
#' @seealso \code{\link{NPEtri}}, \code{\link{NAStri}}, \code{\link{NCStri}},
#' and \code{\link{IarcPEbasic.tri}}
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
#' M<-as.numeric(runif.basic.tri(1,c1,c2)$g)  #try also M<-c(.6,.2)
#'
#' r<-2
#'
#' P1<-as.numeric(runif.basic.tri(1,c1,c2)$g)  #try also P1<-c(.4,.2)
#' NPEbasic.tri(P1,r,c1,c2,M)
#'
#' #or try
#' Rv<-rel.vert.basic.tri(P1,c1,c2,M)$rv
#' NPEbasic.tri(P1,r,c1,c2,M,Rv)
#'
#' P1<-c(1.4,1.2)
#' P2<-c(1.5,1.26)
#' NPEbasic.tri(P1,r,c1,c2,M) #gives an error if M=c(1.3,1.3)
#' #since center is not the circumcenter or not in the interior of the triangle
#' }
#'
#' @export NPEbasic.tri
NPEbasic.tri <- function(p,r,c1,c2,M=c(1,1,1),rv=NULL)
{
  if (!is.point(p) )
  {stop('p must be a numeric 2D point')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  A<-c(0,0); B<-c(1,0); C<-c(c1,c2); Tb<-rbind(A,B,C)
  #standard basic triangle

  CC = circumcenter.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (!in.triangle(p,Tb,boundary=TRUE)$in.tri)
  {reg<-NULL; return(reg); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),
               rel.vert.triCC(p,Tb)$rv,
               rel.vert.tri(p,Tb,M)$rv)
  #vertex region for pt
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}


  x1<-p[1]; y1<-p[2]
  if (rv==1)
  {
    A1<-c(-y1*c1*r/c2+x1*r+y1*r/c2, 0)
    A2<-c(-c1^2*y1*r/c2+c1*x1*r+y1*c1*r/c2, -c1*r*y1+c2*r*x1+r*y1)
    reg<-rbind(A,A1,A2)
  } else if (rv==2)
  {
    B1<-c(-r+c1^2*y1*r/c2-c1*x1*r+c1*r-y1*c1*r/c2+x1*r+1, c1*r*y1-c2*r*x1+c2*r)
    B2<-c(-r-y1*c1*r/c2+x1*r+1, 0)
    reg<-rbind(B,B1,B2)
  } else
  {
    C1<-c(-c1*r+y1*c1*r/c2+c1, -c2*r+r*y1+c2)
    C2<-c(y1*c1*r/c2+r-y1*r/c2-c1*r+c1, -c2*r+r*y1+c2)
    reg<-rbind(C,C1,C2)
  }
  if (abs(area.polygon(reg))>abs(area.polygon(Tb)))
  {reg<-Tb}
  row.names(reg)<-c()
  reg
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))}
#' for points \code{p1} and \code{p2}
#' in the standard basic triangle,
#' that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' and returns 0 otherwise,
#' where \eqn{N_{PE}(x,r)} is the PE proximity region
#' for point \eqn{x} with expansion parameter \eqn{r \ge 1}.
#'
#' PE proximity region is defined
#' with respect to the standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0}
#' and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Vertex regions are based on the center,
#' \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in
#' barycentric coordinates
#' in the interior of the standard basic triangle \eqn{T_b}
#' or based on circumcenter of \eqn{T_b};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \eqn{T_b}.
#' \code{rv} is the index of the vertex region \code{p1} resides,
#' with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct and either of them are
#' outside \eqn{T_b}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations
#' (i.e., it allows loops).
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle.
#' Hence, standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:arc-density-PE;textual}{pcds}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point.
#' The function determines whether \code{p2} is
#' inside the PE proximity region of
#' \code{p1} or not.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}
#' @param c1,c2 Positive real numbers
#'  which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges;
#' \eqn{c_1} must be in \eqn{[0,1/2]},
#' \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard basic triangle
#' or circumcenter of \eqn{T_b}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#' @param rv The index of the vertex region in \eqn{T_b} containing the point,
#' either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))}
#' for points \code{p1} and \code{p2} in the standard basic triangle,
#' that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{IarcPEtri}} and \code{\link{IarcPEstd.tri}}
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
#' r<-2
#'
#' P1<-as.numeric(runif.basic.tri(1,c1,c2)$g)
#' P2<-as.numeric(runif.basic.tri(1,c1,c2)$g)
#' IarcPEbasic.tri(P1,P2,r,c1,c2,M)
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' IarcPEbasic.tri(P1,P2,r,c1,c2,M)
#' IarcPEbasic.tri(P2,P1,r,c1,c2,M)
#'
#' #or try
#' Rv<-rel.vert.basic.tri(P1,c1,c2,M)$rv
#' IarcPEbasic.tri(P1,P2,r,c1,c2,M,Rv)
#' }
#'
#' @export IarcPEbasic.tri
IarcPEbasic.tri <- function(p1,p2,r,c1,c2,M=c(1,1,1),rv=NULL)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  CC = circumcenter.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,Tb,boundary=TRUE)$in.tri ||
      !in.triangle(p2,Tb,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),
               rel.vert.triCC(p1,Tb)$rv,rel.vert.tri(p1,Tb,M)$rv)
  #vertex region for p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  X1<-p1[1]; Y1<-p1[2];
  X2<-p2[1]; Y2<-p2[2];
  arc<-0;
  if (rv==1)
  {
    x1n<-X1*r; y1n<-Y1*r;
    if ( Y2 < paraline(c(x1n,y1n),y2,y3,X2)$y ) {arc <-1}
  } else {
    if (rv==2)
    {
      x1n<-1+(X1-1)*r; y1n<-Y1*r;
      if ( Y2 < paraline(c(x1n,y1n),y1,y3,X2)$y ) {arc <-1}
    } else {
      y1n<-y3[2]+(Y1-y3[2])*r;
      if ( Y2 > y1n ) {arc<-1}
    }}
  arc
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another
#' for Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))}
#' for points \code{p1} and \code{p2}
#' in the standard equilateral triangle,
#' that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' and returns 0 otherwise,
#' where \eqn{N_{PE}(x,r)} is the PE proximity region
#' for point \eqn{x} with expansion parameter \eqn{r \ge 1}.
#'
#' PE proximity region is defined
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_e}.
#' \code{rv} is the index of the vertex region \code{p1} resides,
#' with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct
#' and either of them are outside \eqn{T_e}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations
#' (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:arc-density-CS;textual}{pcds}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether \code{p2} is inside the PE proximity region of
#' \code{p1} or not.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param rv The index of the vertex region in \eqn{T_e}
#' containing the point, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))}
#' for points \code{p1} and \code{p2}
#' in the standard equilateral triangle,
#' that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{IarcPEtri}}, \code{\link{IarcPEbasic.tri}},
#' and \code{\link{IarcCSstd.tri}}
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
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' IarcPEstd.tri(Xp[1,],Xp[3,],r=1.5,M)
#' IarcPEstd.tri(Xp[1,],Xp[3,],r=2,M)
#'
#' #or try
#' Rv<-rel.vert.std.triCM(Xp[1,])$rv
#' IarcPEstd.tri(Xp[1,],Xp[3,],r=2,rv=Rv)
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' r<-2
#' IarcPEstd.tri(P1,P2,r,M)
#' }
#'
#' @export IarcPEstd.tri
IarcPEstd.tri <- function(p1,p2,r,M=c(1,1,1),rv=NULL)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
        3D point for barycentric coordinates')}

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==F)
  {stop('center is not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,Te,boundary=TRUE)$in.tri ||
      !in.triangle(p2,Te,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  {rv<-rel.vert.std.tri(p1,M)$rv #vertex region for p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  X1<-p1[1]; Y1<-p1[2];
  X2<-p2[1]; Y2<-p2[2];
  arc<-0;
  if (rv==1)
  {
    x1n<-X1*r; y1n<-Y1*r;
    if ( Y2 < y1n-3^(1/2)*X2+3^(1/2)*x1n ) {arc <-1}
  } else {
    if (rv==2)
    {
      x1n<-1+(X1-1)*r; y1n<-Y1*r;
      if ( Y2 < y1n+3^(1/2)*X2-3^(1/2)*x1n  ) {arc <-1}
    } else {
      y1n<-C[2]+(Y1-C[2])*r;
      if ( Y2 > y1n ) {arc<-1}
    }}
  arc
} #end of the function
#'

#################################################################

#' @title Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' and quantities related to the triangle - standard equilateral triangle case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}
#' in the standard equilateral triangle.
#' It also provides number of vertices
#' (i.e., number of data points inside the standard equilateral triangle \eqn{T_e})
#' and indices of the data points that reside in \eqn{T_e}.
#'
#' PE proximity region \eqn{N_{PE}(x,r)} is defined
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_e}.
#' For the number of arcs, loops are not allowed so
#' arcs are only possible for points inside \eqn{T_e} for this function.
#'
#' See also (\insertCite{ceyhan:arc-density-PE;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter for PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and quantities related to the standard equilateral triangle}
#' \item{num.arcs}{Number of arcs of the PE-PCD}
#' \item{num.in.tri}{Number of \code{Xp} points
#' in the standard equilateral triangle, \eqn{T_e}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in \eqn{T_e}}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support triangle \eqn{T_e}.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsPEtri}}, \code{\link{num.arcsPE}},
#' and \code{\link{num.arcsCSstd.tri}}
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
#' M<-c(.6,.2)  #try also M<-c(1,1,1)
#'
#' Narcs = num.arcsPEstd.tri(Xp,r=1.25,M)
#' Narcs
#' summary(Narcs)
#' par(pty="s")
#' plot(Narcs,asp=1)
#' }
#'
#' @export num.arcsPEstd.tri
num.arcsPEstd.tri <- function(Xp,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
        3D point for barycentric coordinates')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==F)
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
    {p1<-Xp[i,]
    if (!in.triangle(p1,Te,boundary = TRUE)$in.tri)
    {arcs<-arcs+0
    } else
    {
      ind.in.tri = c(ind.in.tri,i)
      rv<-rel.vert.std.tri(p1,M)$rv

      for (j in ((1:n)[-i]) )
      {p2<-Xp[j,]
      if (!in.triangle(p2,Te,boundary = TRUE)$in.tri)
      {arcs<-arcs+0
      } else
      {
        arcs<-arcs+IarcPEstd.tri(p1,p2,r,M,rv)
      }
      }
    }
    }
  }

  NinTri = length(ind.in.tri)
  desc<-"Number of Arcs of the PE-PCD and the Related Quantities with vertices Xp in the Standard Equilateral Triangle"

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

#' @title Incidence matrix for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - standard equilateral triangle case
#'
#' @description Returns the incidence matrix for the PE-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#'
#' PE proximity region is constructed
#' with respect to the standard equilateral triangle \eqn{T_e} with
#' expansion parameter \eqn{r \ge 1} and vertex regions are based on
#' the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of \eqn{T_e}; default is \eqn{M=(1,1,1)},
#' i.e., the center of mass of \eqn{T_e}.
#' Loops are allowed,
#' so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return Incidence matrix for the PE-PCD with vertices
#' being 2D data set, \code{Xp}
#' in the standard equilateral triangle where PE proximity
#' regions are defined with \code{M}-vertex regions.
#'
#' @seealso \code{\link{inci.matPEtri}}, \code{\link{inci.matPE}},
#' and \code{\link{inci.matCSstd.tri}}
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
#' n<-10
#'
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(runif.std.tri(1)$g)  #try also M<-c(.6,.2)
#'
#' inc.mat<-inci.matPEstd.tri(Xp,r=1.25,M)
#' inc.mat
#' sum(inc.mat)-n
#' num.arcsPEstd.tri(Xp,r=1.25)
#'
#' dom.num.greedy(inc.mat)
#' Idom.num.up.bnd(inc.mat,2) #try also dom.num.exact(inc.mat)
#' }
#'
#' @export inci.matPEstd.tri
inci.matPEstd.tri <- function(Xp,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==F)
  {stop('center is not in the interior of the triangle')}

  n<-nrow(Xp)

  inc.mat<-matrix(0, nrow=n, ncol=n)
  for (i in 1:n)
  {p1<-Xp[i,]
  rv<-rel.vert.std.tri(p1,M)$rv
  for (j in ((1:n)) )
  {p2<-Xp[j,]
  inc.mat[i,j]<-IarcPEstd.tri(p1,p2,r,M,rv)
  }
  }
  inc.mat
} #end of the function
#'

#################################################################

# funsMuVarPE2D
#'
#' @title Returns the mean and (asymptotic) variance of arc density of
#' Proportional Edge Proximity Catch Digraph (PE-PCD)
#' for 2D uniform data in one triangle
#'
#' @description
#' Two functions: \code{muPE2D} and \code{asyvarPE2D}.
#'
#' \code{muPE2D} returns the mean of the (arc) density of PE-PCD
#' and \code{asyvarPE2D} returns the asymptotic variance
#' of the arc density of PE-PCD
#' for 2D uniform data in a triangle.
#'
#' PE proximity regions are defined
#' with expansion parameter \eqn{r \ge 1}
#' with respect to the triangle
#' in which the points reside and
#' vertex regions are based on center of mass, \eqn{CM} of the triangle.
#'
#' See also (\insertCite{ceyhan:arc-density-PE;textual}{pcds}).
#'
#' @param r A positive real number which serves
#' as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#'
#' @return \code{muPE2D} returns the mean
#' and \code{asyvarPE2D} returns the (asymptotic) variance of the
#' arc density of PE-PCD for uniform data in any triangle.
#'
#' @name funsMuVarPE2D
NULL
#'
#' @seealso \code{\link{muCS2D}} and \code{\link{asyvarCS2D}}
#'
#' @rdname funsMuVarPE2D
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for muPE2D
#' muPE2D(1.2)
#'
#' rseq<-seq(1.01,5,by=.05)
#' lrseq<-length(rseq)
#'
#' mu<-vector()
#' for (i in 1:lrseq)
#' {
#'   mu<-c(mu,muPE2D(rseq[i]))
#' }
#'
#' plot(rseq, mu,type="l",xlab="r",ylab=expression(mu(r)),lty=1,
#' xlim=range(rseq),ylim=c(0,1))
#' }
#'
#' @export muPE2D
muPE2D <- function(r)
{
  if (!is.point(r,1) || r<1)
  {stop('The argument must be a scalar greater than 1')}

  mn<-0;
  if (r < 3/2)
  {
    mn<-(37/216)*r^2;
  } else {
    if (r < 2)
    {
      mn<--r^2/8-8/r+9/(2*r^2)+4;
    } else {
      mn<-1-3/(2*r^2);
    }}
  mn
} #end of the function
#'
#' @references
#' \insertAllCited{}
#'
#' @rdname funsMuVarPE2D
#'
#' @examples
#' \dontrun{
#' #Examples for asyvarPE2D
#' asyvarPE2D(1.2)
#'
#' rseq<-seq(1.01,5,by=.05)
#' lrseq<-length(rseq)
#'
#' avar<-vector()
#' for (i in 1:lrseq)
#' {
#'   avar<-c(avar,asyvarPE2D(rseq[i]))
#' }
#'
#' par(mar=c(5,5,4,2))
#' plot(rseq, avar,type="l",xlab="r",
#' ylab=expression(paste(sigma^2,"(r)")),lty=1,xlim=range(rseq))
#' }
#'
#' @export asyvarPE2D
asyvarPE2D <- function(r)
{
  if (!is.point(r,1) || r<1)
  {stop('The argument must be a scalar greater than 1')}

  asyvar<-0;
  if (r < 4/3)
  {
    asyvar<-(3007*r^(10)-13824*r^9+898*r^8+77760*r^7-117953*r^6+48888*r^5-24246*r^4+60480*r^3-38880*r^2+3888)/(58320*r^4);
  } else {
    if (r < 3/2)
    {
      asyvar<-(5467*r^(10)-37800*r^9+61912*r^8+46588*r^6-191520*r^5+13608*r^4+241920*r^3-155520*r^2+15552)/(233280*r^4);
    } else {
      if (r < 2)
      {
        asyvar<--(7*r^(12)-72*r^(11)+312*r^(10)-5332*r^8+15072*r^7+13704*r^6-139264*r^5+273600*r^4-242176*r^3+103232*r^2-27648*r+8640)/(960*r^6);
      } else {
        asyvar<-(15*r^4-11*r^2-48*r+25)/(15*r^6);
      }}}
  asyvar # no need to multiply this by 4
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point or not for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p} is a dominating point
#' of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp}
#' for data in the standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' that is, returns 1 if \code{p} is a dominating point of PE-PCD,
#' and returns 0 otherwise.
#'
#' PE proximity regions are defined
#' with respect to the standard basic triangle \eqn{T_b}.
#' In the standard basic triangle, \eqn{T_b},
#' \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the
#' original triangle. Hence, standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)} in Cartesian
#' coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of a standard basic triangle
#' to the edges on the extension of the lines joining \code{M}
#' to the vertices or based on the circumcenter of \eqn{T_b};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#' Point, \code{p}, is in the vertex region of vertex \code{rv}
#' (default is \code{NULL});
#' vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise.
#'
#' \code{ch.data.pnt} is for checking
#' whether point \code{p} is a data point in \code{Xp} or not
#' (default is \code{FALSE}),
#' so by default this function checks
#' whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point
#' or not of the PE-PCD.
#' @param Xp A set of 2D points
#' which constitutes the vertices of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c1,c2 Positive real numbers
#' which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]},
#' \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard basic triangle \eqn{T_b}
#' or the circumcenter of \eqn{T_b}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#' @param rv Index of the vertex whose region contains point \code{p},
#' \code{rv} takes the vertex labels as \eqn{1,2,3} as
#' in the row order of the vertices in \eqn{T_b}.
#' @param ch.data.pnt A logical argument for checking
#' whether point \code{p} is a data point in \code{Xp} or not
#' (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point,
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{Idom.num1ASbasic.tri}} and \code{\link{Idom.num1AStri}}
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
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.basic.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.basic.tri(1,c1,c2)$g)  #try also M<-c(.6,.3)
#' r<-2
#'
#' P<-c(.4,.2)
#' Idom.num1PEbasic.tri(P,Xp,r,c1,c2,M)
#' Idom.num1PEbasic.tri(Xp[1,],Xp,r,c1,c2,M)
#'
#' Idom.num1PEbasic.tri(c(1,1),Xp,r,c1,c2,M,ch.data.pnt = FALSE)
#' #gives an error message if ch.data.pnt = TRUE since point p=c(1,1) is not a data point in Xp
#'
#' #or try
#' Rv<-rel.vert.basic.tri(Xp[1,],c1,c2,M)$rv
#' Idom.num1PEbasic.tri(Xp[1,],Xp,r,c1,c2,M,Rv)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1PEbasic.tri(Xp[i,],Xp,r,c1,c2,M))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tb)}
#' #need to run this when M is given in barycentric coordinates
#'
#' if (identical(M,circumcenter.tri(Tb)))
#' {
#'   plot(Tb,pch=".",asp=1,xlab="",ylab="",axes=TRUE,
#'   xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#'   polygon(Tb)
#'   points(Xp,pch=1,col=1)
#'   Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2)
#' } else
#' {plot(Tb,pch=".",xlab="",ylab="",axes=TRUE,
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#'   polygon(Tb)
#'   points(Xp,pch=1,col=1)
#'   Ds<-prj.cent2edges.basic.tri(c1,c2,M)}
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(rbind(Xp[ind.gam1,]),pch=4,col=2)
#'
#' txt<-rbind(Tb,M,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,-.02,.03,-.03,.01)
#' yc<-txt[,2]+c(.02,.02,.02,-.02,.02,.02,-.03)
#' txt.str<-c("A","B","C","M","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' Idom.num1PEbasic.tri(c(.2,.1),Xp,r,c1,c2,M,ch.data.pnt=FALSE)
#' #gives an error message if ch.data.pnt=TRUE since point p is not a data point in Xp
#' }
#'
#' @export Idom.num1PEbasic.tri
Idom.num1PEbasic.tri <- function(p,Xp,r,c1,c2,M=c(1,1,1),rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('p is not a data point in Xp')}
  }

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  CC = circumcenter.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (in.triangle(p,Tb)$in.tri==F)
  {dom<-0; return(dom); stop}
  #('point is not inside the triangle')}

  n<-nrow(Xp)
  dom<-1; i<-1;

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rel.vert.triCC(p,Tb)$rv,
               rel.vert.tri(p,Tb,M)$rv)  #vertex region for p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  while (i <= n & dom==1)
  {if (IarcPEbasic.tri(p,Xp[i,],r,c1,c2,M,rv)==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for two points being a dominating set for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard basic triangle case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is
#' a dominating set of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp}
#' in the standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))},
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD,
#' and returns 0 otherwise.
#'
#' PE proximity regions are defined with respect to \eqn{T_b}.
#' In the standard basic triangle, \eqn{T_b},
#' \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0}
#' and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the
#' original triangle.
#' Hence, standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)} in Cartesian
#' coordinates or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of a standard basic triangle \eqn{T_b};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#' Point, \code{p1},
#' is in the vertex region of vertex \code{rv1}
#' (default is \code{NULL});
#' and point, \code{p2},
#' is in the vertex region of vertex \code{rv2}
#' (default is \code{NULL});
#' vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise.
#'
#' \code{ch.data.pnts} is for checking
#' whether points \code{p1} and \code{p2} are both data points
#' in \code{Xp} or not (default is \code{FALSE}),
#' so by default this function checks
#' whether the points \code{p1} and \code{p2} would constitute
#' a dominating set
#' if they both were actually in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param p1,p2 Two 2D points to be tested for
#' constituting a dominating set of the PE-PCD.
#' @param Xp A set of 2D points
#' which constitutes the vertices of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param c1,c2 Positive real numbers
#' which constitute the vertex of the standard basic triangle.
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]},
#' \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard basic triangle \eqn{T_b}
#' or the circumcenter of \eqn{T_b}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_b}.
#' @param rv1,rv2 The indices of the vertices
#' whose regions contains \code{p1} and \code{p2}, respectively.
#' They take the vertex labels as \eqn{1,2,3} as
#' in the row order of the vertices in \eqn{T_b}
#' (default is \code{NULL} for both).
#' @param ch.data.pnts A logical argument for
#' checking whether points \code{p1} and \code{p2}
#'  are data points in \code{Xp} or not
#' (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD,
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{Idom.num2PEtri}}, \code{\link{Idom.num2ASbasic.tri}},
#' and \code{\link{Idom.num2AStri}}
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
#' n<-10  #try also n<-20
#'
#' set.seed(1)
#' Xp<-runif.basic.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.basic.tri(1,c1,c2)$g)  #try also M<-c(.6,.3)
#'
#' r<-2
#'
#' Idom.num2PEbasic.tri(Xp[1,],Xp[2,],Xp,r,c1,c2,M)
#'
#' Idom.num2PEbasic.tri(c(1,2),c(1,3),rbind(c(1,2),c(1,3)),r,c1,c2,M)
#' Idom.num2PEbasic.tri(c(1,2),c(1,3),rbind(c(1,2),c(1,3)),r,c1,c2,M,
#' ch.data.pnts = TRUE)
#'
#' ind.gam2<-vector()
#' for (i in 1:(n-1))
#'   for (j in (i+1):n)
#'   {if (Idom.num2PEbasic.tri(Xp[i,],Xp[j,],Xp,r,c1,c2,M)==1)
#'    ind.gam2<-rbind(ind.gam2,c(i,j))}
#' ind.gam2
#'
#' #or try
#' rv1<-rel.vert.basic.tri(Xp[1,],c1,c2,M)$rv;
#' rv2<-rel.vert.basic.tri(Xp[2,],c1,c2,M)$rv;
#' Idom.num2PEbasic.tri(Xp[1,],Xp[2,],Xp,r,c1,c2,M,rv1,rv2)
#'
#' #or try
#' rv1<-rel.vert.basic.tri(Xp[1,],c1,c2,M)$rv;
#' Idom.num2PEbasic.tri(Xp[1,],Xp[2,],Xp,r,c1,c2,M,rv1)
#'
#' #or try
#' rv2<-rel.vert.basic.tri(Xp[2,],c1,c2,M)$rv;
#' Idom.num2PEbasic.tri(Xp[1,],Xp[2,],Xp,r,c1,c2,M,rv2=rv2)
#'
#' Idom.num2PEbasic.tri(c(1,2),Xp[2,],Xp,r,c1,c2,M,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE since not both points are data points in Xp
#' }
#'
#' @export Idom.num2PEbasic.tri
Idom.num2PEbasic.tri <- function(p1,p2,Xp,r,c1,c2,M=c(1,1,1),rv1=NULL,rv2=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points are data points in Xp')}
  }

  if (isTRUE(all.equal(p1,p2)))
  {dom<-0; return(dom); stop}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  CC = circumcenter.tri(Tb)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,Tb,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (is.null(rv1))
  {rv1<-ifelse(isTRUE(all.equal(M,CC)),rel.vert.triCC(p1,Tb)$rv,
               rel.vert.tri(p1,Tb,M)$rv)  #vertex region for point p1
  }
  if (is.null(rv2))
  {rv2<-ifelse(isTRUE(all.equal(M,CC)),rel.vert.triCC(p2,Tb)$rv,
               rel.vert.tri(p2,Tb,M)$rv)  #vertex region for point p2
  }

  Xp<-matrix(Xp,ncol=2)
  n<-nrow(Xp)
  dom<-1; i<-1;

  while (i <= n & dom==1)
  {if (max(IarcPEbasic.tri(p1,Xp[i,],r,c1,c2,M,rv1),
           IarcPEbasic.tri(p2,Xp[i,],r,c1,c2,M,rv2))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The vertices of the Proportional Edge (PE) Proximity Region
#' in a general triangle
#'
#' @description Returns the vertices of the PE proximity region
#' (which is itself a triangle) for a point in the
#' triangle \code{tri}\eqn{=T(A,B,C)=}\code{(rv=1,rv=2,rv=3)}.
#'
#' PE proximity region is defined with respect to the triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#'
#' Vertex regions are labeled as \eqn{1,2,3}
#' rowwise for the vertices
#' of the triangle \code{tri}.
#' \code{rv} is the index of the vertex region \code{p} resides,
#' with default=\code{NULL}.
#' If \code{p} is outside of \code{tri},
#' it returns \code{NULL} for the proximity region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param p A 2D point whose PE proximity region is to be computed.
#' @param r A positive real number which serves
#' as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param rv Index of the \code{M}-vertex region
#' containing the point \code{p}, either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return Vertices of the triangular region
#' which constitutes the PE proximity region with expansion parameter
#' \code{r} and center \code{M} for a point \code{p}
#'
#' @seealso \code{\link{NPEbasic.tri}}, \code{\link{NAStri}},
#' \code{\link{NCStri}}, and \code{\link{IarcPEtri}}
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
#' r<-1.5
#'
#' n<-3
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' NPEtri(Xp[3,],Tr,r,M)
#'
#' P1<-as.numeric(runif.tri(1,Tr)$g)  #try also P1<-c(.4,.2)
#' NPEtri(P1,Tr,r,M)
#'
#' M<-c(1.3,1.3)
#' r<-2
#'
#' P1<-c(1.4,1.2)
#' P2<-c(1.5,1.26)
#' NPEtri(P1,Tr,r,M)
#' NPEtri(P2,Tr,r,M)
#'
#' #or try
#' Rv<-rel.vert.tri(P1,Tr,M)$rv
#' NPEtri(P1,Tr,r,M,Rv)
#' }
#'
#' @export NPEtri
NPEtri <- function(p,tri,r,M=c(1,1,1),rv=NULL)
{
  if (!is.point(p) )
  {stop('p must be a numeric 2D point')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (!in.triangle(p,tri,boundary=TRUE)$in.tri)
  {reg<-NULL; return(reg); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rel.vert.triCC(p,tri)$rv,
               rel.vert.tri(p,tri,M)$rv)  #vertex region for p
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  A<-tri[1,];  B<-tri[2,];  C<-tri[3,];
  if (rv==1)
  {
    d1<-dist.point2line(p,B,C)$dis
    d2<-dist.point2line(A,B,C)$dis
    sr<-d1/d2
    P1<-B+sr*(A-B); P2<-C+sr*(A-C)
    A1<-A+r*(P1-A)
    A2<-A+r*(P2-A)
    reg<-rbind(A,A1,A2)
  } else if (rv==2)
  {
    d1<-dist.point2line(p,A,C)$dis
    d2<-dist.point2line(B,A,C)$dis
    sr<-d1/d2
    P1<-A+sr*(B-A); P2<-C+sr*(B-C)
    B1<-B+r*(P1-B)
    B2<-B+r*(P2-B)
    reg<-rbind(B,B1,B2)
  } else
  {
    d1<-dist.point2line(p,A,B)$dis
    d2<-dist.point2line(C,A,B)$dis
    sr<-d1/d2
    P1<-A+sr*(C-A); P2<-B+sr*(C-B)
    C1<-C+r*(P1-C)
    C2<-C+r*(P2-C)
    reg<-rbind(C,C1,C2)
  }
  if (abs(area.polygon(reg))>abs(area.polygon(tri)))
  {reg<-tri}
  row.names(reg)<-c()
  reg
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point to another
#' for Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description Returns \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))}
#' for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' and returns 0 otherwise,
#' where \eqn{N_{PE}(x,r)} is the PE proximity region for point \eqn{x}
#' with the expansion parameter \eqn{r \ge 1}.
#'
#' PE proximity region is constructed
#' with respect to the triangle \code{tri} and
#' vertex regions are based on the center, \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' \code{rv} is the index of the vertex region \code{p1} resides,
#' with default=\code{NULL}.
#'
#' If \code{p1} and \code{p2} are distinct
#' and either of them are outside \code{tri}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations
#' (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point.
#' The function determines whether \code{p2} is
#' inside the PE proximity region of
#' \code{p1} or not.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param rv Index of the \code{M}-vertex region containing the point,
#' either \code{1,2,3} or \code{NULL}
#' (default is \code{NULL}).
#'
#' @return \eqn{I(}\code{p2} is in \eqn{N_{PE}(p1,r))}
#' for points \code{p1} and \code{p2},
#' that is, returns 1 if \code{p2} is in \eqn{N_{PE}(p1,r)},
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{IarcPEbasic.tri}}, \code{\link{IarcPEstd.tri}},
#' \code{\link{IarcAStri}}, and \code{\link{IarcCStri}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0);
#'
#' r<-1.5
#'
#' n<-3
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' IarcPEtri(Xp[1,],Xp[2,],Tr,r,M)
#'
#' P1<-as.numeric(runif.tri(1,Tr)$g)
#' P2<-as.numeric(runif.tri(1,Tr)$g)
#' IarcPEtri(P1,P2,Tr,r,M)
#'
#' P1<-c(.4,.2)
#' P2<-c(1.8,.5)
#' IarcPEtri(P1,P2,Tr,r,M)
#' IarcPEtri(P2,P1,Tr,r,M)
#'
#' M<-c(1.3,1.3)
#' r<-2
#'
#' #or try
#' Rv<-rel.vert.tri(P1,Tr,M)$rv
#' IarcPEtri(P1,P2,Tr,r,M,Rv)
#' }
#'
#' @export IarcPEtri
IarcPEtri <- function(p1,p2,tri,r,M=c(1,1,1),rv=NULL)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {arc<-1; return(arc); stop}

  if (!in.triangle(p1,tri,boundary=TRUE)$in.tri ||
      !in.triangle(p2,tri,boundary=TRUE)$in.tri)
  {arc<-0; return(arc); stop}

  if (is.null(rv))
  { rv<-ifelse(isTRUE(all.equal(M,CC)),rel.vert.triCC(p1,tri)$rv,
               rel.vert.tri(p1,tri,M)$rv)  #vertex region for p1
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  pr<-NPEtri(p1,tri,r,M,rv)  #proximity region
  arc<-ifelse(area.polygon(pr)==0,0,
              sum(in.triangle(p2,pr,boundary=TRUE)$in.tri))

  arc
} #end of the function
#'

#################################################################

#' @title Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' and quantities related to the triangle - one triangle case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' PE proximity region \eqn{N_{PE}(x,r)} is defined
#' with respect to the triangle, \code{tri}
#' with expansion parameter \eqn{r \ge 1} and vertex regions are
#' based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri} or
#' based on circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' For the number of arcs, loops are not
#' allowed so arcs are only possible for points
#' inside the triangle \code{tri} for this function.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and quantities related to the triangle}
#' \item{num.arcs}{Number of arcs of the PE-PCD}
#' \item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in the triangle}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the support triangle.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsPEstd.tri}}, \code{\link{num.arcsPE}},
#' \code{\link{num.arcsCStri}}, and \code{\link{num.arcsAStri}}
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
#' Narcs = num.arcsPEtri(Xp,Tr,r=1.25,M)
#' Narcs
#' summary(Narcs)
#' plot(Narcs)
#' }
#'
#' @export num.arcsPEtri
num.arcsPEtri <- function(Xp,tri,r,M=c(1,1,1))
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  n<-nrow(Xp)
  arcs<-0
  ind.in.tri = NULL
  if (n<=0)
  {
    arcs<-0
  } else
  {
    for (i in 1:n)
    {
      if (in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri)
      {  vert<-ifelse(isTRUE(all.equal(M,CC)),
                      rel.vert.triCC(Xp[i,],tri)$rv,
                      rel.vert.tri(Xp[i,],tri,M)$rv)
      ind.in.tri = c(ind.in.tri,i)
      for (j in (1:n)[-i])  #to avoid loops
      {
        arcs<-arcs+IarcPEtri(Xp[i,],Xp[j,],tri,r,M,rv=vert)
      }
      }
    }
  }

  NinTri = length(ind.in.tri)
  desc<-"Number of Arcs of the PE-PCD with vertices Xp and Quantities Related to the Support Triangle"

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

#' @title Arc density of Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - one triangle case
#'
#' @description Returns the arc density of PE-PCD
#' whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the triangle \code{tri}.
#'
#' PE proximity regions is defined with respect to \code{tri} with
#' expansion parameter \eqn{r \ge 1} and vertex regions are
#' based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri} or based on
#' circumcenter of \code{tri}; default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' The function also provides arc density standardized
#' by the mean and asymptotic variance of the arc density
#' of PE-PCD for uniform data in the triangle \code{tri}
#' only when \code{M} is the center of mass.
#' For the number of arcs, loops are not allowed.
#'
#' \code{tri.cor} is a logical argument for triangle correction
#' (default is \code{TRUE}),
#' if \code{TRUE}, only the points
#' inside the triangle are considered
#' (i.e., digraph induced by these vertices are considered) in computing
#' the arc density, otherwise all points are considered
#' (for the number of vertices in the denominator of arc density).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param tri.cor A logical argument
#' for computing the arc density for only the points inside the triangle,
#' \code{tri}.
#' (default is \code{tri.cor=FALSE}), i.e.,
#' if \code{tri.cor=TRUE} only the induced digraph with the vertices
#' inside \code{tri} are considered in the
#' computation of arc density.
#'
#' @return A \code{list} with the elements
#' \item{arc.dens}{Arc density of PE-PCD
#' whose vertices are the 2D numerical data set, \code{Xp};
#' PE proximity regions are defined
#' with respect to the triangle \code{tri} and \code{M}-vertex regions}
#' \item{std.arc.dens}{Arc density standardized
#' by the mean and asymptotic variance of the arc
#' density of PE-PCD for uniform data in the triangle \code{tri}.
#' This will only be returned, if \code{M} is the center of mass.}
#'
#' @seealso \code{\link{ASarc.dens.tri}}, \code{\link{CSarc.dens.tri}},
#' and \code{\link{num.arcsPEtri}}
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
#' num.arcsPEtri(Xp,Tr,r=1.5,M)
#' PEarc.dens.tri(Xp,Tr,r=1.5,M)
#' PEarc.dens.tri(Xp,Tr,r=1.5,M,tri.cor = TRUE)
#' }
#'
#' @export PEarc.dens.tri
PEarc.dens.tri <- function(Xp,tri,r,M=c(1,1,1),tri.cor=FALSE)
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  nx<-nrow(Xp)

  narcs<-num.arcsPEtri(Xp,tri,r,M)$num.arcs
  mean.rho<-muPE2D(r)
  var.rho<-asyvarPE2D(r)

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
    {stop('There are not enough Xp points in the triangle, tri,
          to compute the arc density!')}
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

#################################################################

#' @title Number of arcs of Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' and related quantities of the induced subdigraphs for points in the Delaunay triangles -
#' multiple triangle case
#'
#' @description
#' An object of class \code{"NumArcs"}.
#' Returns the number of arcs and various other quantities related to the Delaunay triangles
#' for Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' PE proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions in each triangle
#' is based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each
#' Delaunay triangle or based on circumcenter of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
#' Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points). For the number of arcs,
#' loops are not allowed so arcs are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE;textual}{pcds})
#' for more on PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of arcs
#' and related quantities for the induced subdigraphs in the Delaunay triangles}
#' \item{num.arcs}{Total number of arcs in all triangles,
#' i.e., the number of arcs for the entire PE-PCD}
#' \item{num.in.conhull}{Number of \code{Xp} points
#' in the convex hull of \code{Yp} points}
#' \item{num.in.tris}{The vector of number of \code{Xp} points
#' in the Delaunay triangles based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the areas of
#' Delaunay triangles based on \code{Yp} points}
#' \item{tri.num.arcs}{The \code{vector} of the number of arcs
#' of the component of the PE-PCD in the
#' Delaunay triangles based on \code{Yp} points}
#' \item{del.tri.ind}{A matrix of indices of vertices of
#' the Delaunay triangles based on \code{Yp} points,
#' each column corresponds to the vector of
#' indices of the vertices of one triangle.}
#' \item{data.tri.ind}{A \code{vector} of indices of vertices of
#' the Delaunay triangles in which data points reside,
#' i.e., column number of \code{del.tri.ind} for each \code{Xp} point.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed, here, tessellation
#' is the Delaunay triangulation based on \code{Yp} points.}
#' \item{vertices}{Vertices of the digraph, \code{Xp}.}
#'
#' @seealso \code{\link{num.arcsPEtri}}, \code{\link{num.arcsPEstd.tri}},
#' \code{\link{num.arcsCS}}, and \code{\link{num.arcsAS}}
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
#' Xp<-cbind(runif(nx),runif(nx))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' Narcs = num.arcsPE(Xp,Yp,r=1.25,M)
#' Narcs
#' summary(Narcs)
#' plot(Narcs)
#' }
#'
#' @export num.arcsPE
num.arcsPE <- function(Xp,Yp,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if ((!is.point(M,3) && M!="CC") || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates or
        "CC" for circumcenter')}

  nx<-nrow(Xp); ny<-nrow(Yp)

  #Delaunay triangulation of Yp points
  Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
  Ytri<-matrix(interp::triangles(Ytrimesh)[,1:3],ncol=3);
  #the indices of the vertices of the Delaunay triangles (row-wise)
  ndt<-nrow(Ytri)  #number of Delaunay triangles

  inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
  NinCH<-sum(inCH)  #number of points in the convex hull

  Wvec=vector()
  for (i in 1:ndt)
  {
    ifelse(ndt==1,Tri<-Yp[Ytri,],Tri<-Yp[Ytri[i,],])
    #vertices of ith triangle
    Wvec<-c(Wvec,area.polygon(Tri))
  }

  if (ny==3)
  { tri<-as.basic.tri(Yp)$tri
  NumArcs = num.arcsPEtri(Xp,tri,r,M)
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

  desc<-"Number of Arcs of the PE-PCD with vertices Xp and Related Quantities for the Induced Subdigraph for the Points in the Delaunay Triangle"
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

        ifelse(identical(M,"CC"),cent<-circumcenter.tri(tri),cent<-M)
        num.arcs<-num.arcsPEtri(Xpi,tri,r,cent)$num.arcs
        arcs<-c(arcs,num.arcs)  #number of arcs in all triangles as a vector

      }

      Tot.Arcs<-sum(arcs)  #the total number of arcs in all triangles
    }

    desc<-"Number of Arcs of the PE-PCD with vertices Xp and Related Quantities for the Induced Subdigraphs for the Points in the Delaunay Triangles"

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

#' @title A test of segregation/association
#' based on arc density of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function
#' which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster
#' away from \code{Yp} points) and association
#' (where \code{Xp} points cluster around
#' \code{Yp} points) based on the normal approximation
#' of the arc density of the PE-PCD for uniform 2D data.
#'
#' The function yields the test statistic,
#' \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval, estimate
#' and null value for the parameter of interest
#' (which is the arc density),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points, arc density
#' of PE-PCD whose vertices are \code{Xp} points equals
#' to its expected value under the uniform distribution and
#' \code{alternative} could be two-sided, or left-sided
#' (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated
#' around the centers of the triangles,
#' or segregation).
#'
#' PE proximity region is constructed
#' with the expansion parameter \eqn{r \ge 1} and \eqn{CM}-vertex regions
#' (i.e., the test is not available for a general center \eqn{M}
#' at this version of the function).
#'
#' **Caveat:** This test is currently a conditional test,
#' where \code{Xp} points are assumed to be random,
#' while \code{Yp} points are
#' assumed to be fixed (i.e., the test is conditional on \code{Yp} points).
#' Furthermore,
#' the test is a large sample test when \code{Xp} points
#' are substantially larger than \code{Yp} points,
#' say at least 5 times more.
#' This test is more appropriate when supports of \code{Xp}
#' and \code{Yp} have a substantial overlap.
#' Currently, the \code{Xp} points
#' outside the convex hull of \code{Yp} points
#' are handled with a convex hull correction factor
#' (see the description below and the function code.)
#' However, in the special case of no \code{Xp} points
#' in the convex hull of \code{Yp} points,
#' arc density is taken to be 1,
#' as this is clearly a case of segregation.
#' Removing the conditioning and extending it to
#' the case of non-concurring supports is
#' an ongoing topic of research of the author of the package.
#'
#' \code{ch.cor} is for convex hull correction
#' (default is \code{"no convex hull correction"}, i.e., \code{ch.cor=FALSE})
#' which is recommended
#' when both \code{Xp} and \code{Yp} have the same rectangular support.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE;textual}{pcds})
#' for more on the test based on the arc density of PE-PCDs.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param ch.cor A logical argument for convex hull correction,
#' default \code{ch.cor=FALSE},
#' recommended when both \code{Xp} and \code{Yp}
#' have the same rectangular support.
#' @param alternative Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval,
#' default is \code{0.95}, for the arc density of PE-PCD based on
#' the 2D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test
#' for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval for the arc density
#' at the given confidence level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate of the parameter, i.e., arc density}
#' \item{null.value}{Hypothesized value for the parameter,
#' i.e., the null arc density, which is usually the
#' mean arc density under uniform distribution.}
#' \item{alternative}{Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{CSarc.dens.test}} and \code{\link{PEarc.dens.test1D}}
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
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' plotDelaunay.tri(Xp,Yp,xlab="",ylab="")
#'
#' PEarc.dens.test(Xp,Yp,r=1.25)
#' PEarc.dens.test(Xp,Yp,r=1.25,ch=TRUE)
#' #since Y points are not uniform, convex hull correction is invalid here
#' }
#'
#' @export PEarc.dens.test
PEarc.dens.test <- function(Xp,Yp,r,ch.cor=FALSE,alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one of \"greater\", \"less\", \"two.sided\"")

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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) ||
        conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  Arcs<-num.arcsPE(Xp,Yp,r,M=c(1,1,1))
  #use the default, i.e., CM for the center M
  NinCH<-Arcs$num.in.con

  num.arcs<-Arcs$num.arcs #total number of arcs in the PE-PCD
  num.arcs.tris = Arcs$tri.num.arcs
  #vector of number of arcs in the Delaunay triangles
  num.dat.tris = Arcs$num.in.tris
  #vector of number of data points in the Delaunay triangles
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
      npe = NPEtri(Xpi,tri,r)
      num.arcs = num.arcs+area.polygon(npe)/Wvec[i]
    }
  }
  asy.mean0<-muPE2D(r)  #asy mean value for the r value
  asy.mean<-asy.mean0*sum(LW^2)

  asy.var0<-asyvarPE2D(r)  #asy variance value for the r value
  asy.var<-asy.var0*sum(LW^3)+4*asy.mean0^2*(sum(LW^3)-(sum(LW^2))^2)

  n<-nrow(Xp)  #number of X points
  if (NinCH == 0) {
  warning('There is no Xp point in the convex hull of Yp points to compute arc density,
           but as this is clearly a segregation pattern, arc density is taken to be 1!')
    arc.dens=1
    TS0<-sqrt(n)*(arc.dens-asy.mean)/sqrt(asy.var)  #standardized test stat
  } else
  {  arc.dens<-num.arcs/(NinCH*(NinCH-1))
  TS0<-sqrt(NinCH)*(arc.dens-asy.mean)/sqrt(asy.var)
  #standardized test stat}  #arc density
  }
  estimate1<-arc.dens
  estimate2<-asy.mean

  method <-c("Large Sample z-Test Based on Arc Density of PE-PCD for Testing Uniformity of 2D Data ---")
  if (ch.cor==FALSE)
  {
    TS<-TS0
    method <-c(method, " without Convex Hull Correction")
  }
  else
  {
    m<-nrow(Yp)  #number of Y points
    NoutCH<-n-NinCH #number of points outside of the convex hull

    prop.out<-NoutCH/n #observed proportion of points outside convex hull
    exp.prop.out<-1.7932/m+1.2229/sqrt(m)
    #expected proportion of points outside convex hull

    TS<-TS0+abs(TS0)*sign(prop.out-exp.prop.out)*(prop.out-exp.prop.out)^2
    method <-c(method, " with Convex Hull Correction")
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

#' @title Incidence matrix for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - one triangle case
#'
#' @description Returns the incidence matrix for the PE-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)}.
#'
#' PE proximity regions are constructed with respect to triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#'
#' @return Incidence matrix for the PE-PCD
#' with vertices being 2D data set, \code{Xp},
#' in the triangle \code{tri} with vertex regions based on center \code{M}
#'
#' @seealso \code{\link{inci.matPE}}, \code{\link{inci.matCStri}},
#'  and  \code{\link{inci.matAStri}}
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
#' IM<-inci.matPEtri(Xp,Tr,r=1.25,M)
#'
#' IM
#' dom.num.greedy(IM) #try also dom.num.exact(IM)
#' Idom.num.up.bnd(IM,3)
#' }
#'
#' @export inci.matPEtri
inci.matPEtri <- function(Xp,tri,r,M=c(1,1,1))
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  n<-nrow(Xp)

  inc.mat<-matrix(0, nrow=n, ncol=n)
  if (n>1)
  {
    for (i in 1:n)
    {p1<-Xp[i,]
    rv<-ifelse(isTRUE(all.equal(M,CC)),
               rel.vert.triCC(p1,tri)$rv,
               rel.vert.tri(p1,tri,M)$rv)

    for (j in ((1:n)) )
    {p2<-Xp[j,]
    inc.mat[i,j]<-IarcPEtri(p1,p2,tri,r,M,rv=rv)
    }
    }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The indicator for a point being a dominating point for
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description Returns \eqn{I(}\code{p} is
#' a dominating point of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp}
#' in the triangle \code{tri}, that is,
#' returns 1 if \code{p} is a dominating point of PE-PCD,
#' and returns 0 otherwise.
#'
#' Point, \code{p}, is in the vertex region of vertex \code{rv}
#' (default is \code{NULL}); vertices are labeled as \eqn{1,2,3}
#' in the order they are stacked row-wise in \code{tri}.
#'
#' PE proximity region is constructed
#' with respect to the triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#'
#' \code{ch.data.pnt} is for checking
#' whether point \code{p} is a data point in \code{Xp}
#' or not (default is \code{FALSE}),
#' so by default this function checks
#' whether the point \code{p} would be a dominating point
#' if it actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p A 2D point that is to be tested for being a dominating point
#' or not of the PE-PCD.
#' @param Xp A set of 2D points
#' which constitutes the vertices of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param rv Index of the vertex whose region contains point \code{p},
#' \code{rv} takes the vertex labels as \eqn{1,2,3} as
#' in the row order of the vertices in \code{tri}.
#' @param ch.data.pnt A logical argument for checking
#' whether point \code{p} is a data point
#' in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\code{p} is a dominating point of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \code{p} is a dominating point,
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{Idom.num1PEbasic.tri}} and \code{\link{Idom.num1AStri}}
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
#' r<-1.5  #try also r<-2
#'
#' Idom.num1PEtri(Xp[1,],Xp,Tr,r,M)
#' Idom.num1PEtri(c(1,2),c(1,2),Tr,r,M)
#' Idom.num1PEtri(c(1,2),c(1,2),Tr,r,M,ch.data.pnt = TRUE)
#'
#' gam.vec<-vector()
#' for (i in 1:n)
#' {gam.vec<-c(gam.vec,Idom.num1PEtri(Xp[i,],Xp,Tr,r,M))}
#'
#' ind.gam1<-which(gam.vec==1)
#' ind.gam1
#'
#' #or try
#' Rv<-rel.vert.tri(Xp[1,],Tr,M)$rv
#' Idom.num1PEtri(Xp[1,],Xp,Tr,r,M,Rv)
#'
#' Ds<-prj.cent2edges(Tr,M)
#'
#' if (dimension(M)==3) {M<-bary2cart(M,Tr)}
#' #need to run this when M is given in barycentric coordinates
#'
#' Xlim<-range(Tr[,1],Xp[,1],M[1])
#' Ylim<-range(Tr[,2],Xp[,2],M[2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",axes=TRUE,
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=1,col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(rbind(Xp[ind.gam1,]),pch=4,col=2)
#' #rbind is to insert the points correctly if there is only one dominating point
#'
#' txt<-rbind(Tr,M,Ds)
#' xc<-txt[,1]+c(-.02,.03,.02,-.02,.04,-.03,.0)
#' yc<-txt[,2]+c(.02,.02,.05,-.03,.04,.06,-.07)
#' txt.str<-c("A","B","C","M","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' P<-c(1.4,1)
#' Idom.num1PEtri(P,P,Tr,r,M)
#' Idom.num1PEtri(Xp[1,],Xp,Tr,r,M)
#'
#' Idom.num1PEtri(c(1,2),Xp,Tr,r,M,ch.data.pnt = FALSE)
#' #gives an error message if ch.data.pnt = TRUE since p is not a data point
#' }
#'
#' @export Idom.num1PEtri
Idom.num1PEtri <- function(p,Xp,tri,r,M=c(1,1,1),rv=NULL,ch.data.pnt=FALSE)
{
  if (!is.point(p))
  {stop('p must be a numeric 2D point')}

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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (ch.data.pnt==TRUE)
  {
    if (!is.in.data(p,Xp))
    {stop('p is not a data point in Xp')}
  }

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (in.triangle(p,tri,boundary=TRUE)$in.tri==FALSE)
  {dom<-0; return(dom); stop}

  if (is.null(rv))
  {rv<-ifelse(isTRUE(all.equal(M,CC)),
              rel.vert.triCC(p,tri)$rv,
              rel.vert.tri(p,tri,M)$rv)
  } else
  {  if (!is.numeric(rv) || sum(rv==c(1,2,3))!=1)
  {stop('vertex index, rv, must be 1, 2 or 3')}}

  n<-nrow(Xp)
  if (n>=1)
  {
    dom<-1; i<-1;
    while (i <= n & dom==1)
    {
      if (IarcPEtri(p,Xp[i,],tri,r,M,rv)==0)
        dom<-0;
      i<-i+1;
    }
  } else
  {
    dom<-0
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for two points constituting a dominating set for
#'  Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#'  one triangle case
#'
#' @description Returns \eqn{I(}\{\code{p1,p2}\} is
#' a dominating set of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD,
#' and returns 0 otherwise.
#'
#' Point, \code{p1}, is in the region of vertex \code{rv1}
#' (default is \code{NULL}) and point, \code{p2}, is
#' in the region of vertex \code{rv2}
#' (default is \code{NULL}); vertices (and hence \code{rv1} and \code{rv2})
#' are labeled as \eqn{1,2,3} in the order they are stacked
#' row-wise in \code{tri}.
#'
#' PE proximity regions are defined
#' with respect to the triangle \code{tri} and vertex regions
#' are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#'
#' \code{ch.data.pnts} is
#' for checking whether points \code{p1} and \code{p2} are data points
#' in \code{Xp} or not
#' (default is \code{FALSE}),
#' so by default this function checks
#' whether the points \code{p1} and \code{p2} would be a
#' dominating set if they actually were in the data set.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param p1,p2 Two 2D points to be tested for constituting
#' a dominating set of the PE-PCD.
#' @param Xp A set of 2D points
#' which constitutes the vertices of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' @param rv1,rv2 The indices of the vertices
#'  whose regions contains \code{p1} and \code{p2}, respectively.
#' They take the vertex labels as \eqn{1,2,3}
#' as in the row order of the vertices in \code{tri}
#' (default is \code{NULL} for both).
#' @param ch.data.pnts A logical argument for
#' checking whether points \code{p1} and \code{p2} are
#' data points in \code{Xp} or not (default is \code{FALSE}).
#'
#' @return \eqn{I(}\{\code{p1,p2}\} is a dominating set of the PE-PCD\eqn{)}
#' where the vertices of the PE-PCD are the 2D data set \code{Xp},
#' that is, returns 1 if \{\code{p1,p2}\} is a dominating set of PE-PCD,
#' and returns 0 otherwise.
#'
#' @seealso \code{\link{Idom.num2PEbasic.tri}}, \code{\link{Idom.num2AStri}},
#' and \code{\link{Idom.num2PEtetra}}
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
#' r<-1.5  #try also r<-2
#'
#' Idom.num2PEtri(Xp[1,],Xp[2,],Xp,Tr,r,M)
#'
#' ind.gam2<-vector()
#' for (i in 1:(n-1))
#'   for (j in (i+1):n)
#'   {if (Idom.num2PEtri(Xp[i,],Xp[j,],Xp,Tr,r,M)==1)
#'    ind.gam2<-rbind(ind.gam2,c(i,j))}
#' ind.gam2
#'
#' #or try
#' rv1<-rel.vert.tri(Xp[1,],Tr,M)$rv;
#' rv2<-rel.vert.tri(Xp[2,],Tr,M)$rv
#' Idom.num2PEtri(Xp[1,],Xp[2,],Xp,Tr,r,M,rv1,rv2)
#'
#' Idom.num2PEtri(Xp[1,],c(1,2),Xp,Tr,r,M,ch.data.pnts = FALSE)
#' #gives an error message if ch.data.pnts = TRUE
#' #since not both points, p1 and p2, are data points in Xp
#' }
#'
#' @export Idom.num2PEtri
Idom.num2PEtri <- function(p1,p2,Xp,tri,r,M=c(1,1,1),rv1=NULL,rv2=NULL,ch.data.pnts=FALSE)
{
  if (!is.point(p1) || !is.point(p2))
  {stop('p1 and p2 must be numeric 2D points')}

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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (ch.data.pnts==TRUE)
  {
    if (!is.in.data(p1,Xp) || !is.in.data(p2,Xp))
    {stop('not both points, p1 and p2, are data points in Xp')}
  }

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  if (isTRUE(all.equal(p1,p2)))
  {dom<-0; return(dom); stop}

  if (is.null(rv1))
  {rv1<-ifelse(isTRUE(all.equal(M,CC)),
               rel.vert.triCC(p1,tri)$rv,
               rel.vert.tri(p1,tri,M)$rv)  #vertex region for point p1
  }
  if (is.null(rv2))
  {rv2<-ifelse(isTRUE(all.equal(M,CC)),
               rel.vert.triCC(p2,tri)$rv,
               rel.vert.tri(p2,tri,M)$rv)  #vertex region for point p2
  }

  n<-nrow(Xp)
  dom<-1; i<-1;
  while (i <= n & dom==1)
  {if (max(IarcPEtri(p1,Xp[i,],tri,r,M,rv1),
           IarcPEtri(p2,Xp[i,],tri,r,M,rv2))==0)
    dom<-0;
  i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The domination number of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) - one triangle case
#'
#' @description Returns the domination number of PE-PCD
#' whose vertices are the data points in \code{Xp}.
#'
#' PE proximity region is defined
#' with respect to the triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1} and
#' vertex regions are constructed with center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the digraph.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{(1,1,1)}, i.e., the center of mass.
#'
#' @return A \code{list} with two elements
#' \item{dom.num}{Domination number of PE-PCD with vertex set = \code{Xp}
#' and expansion parameter \eqn{r \ge 1} and center \code{M}}
#' \item{mds}{A minimum dominating set of PE-PCD with vertex set = \code{Xp}
#' and expansion parameter \eqn{r \ge 1} and center \code{M}}
#' \item{ind.mds}{Indices of the minimum dominating set \code{mds}}
#'
#' @seealso \code{\link{PEdom.num.nondeg}}, \code{\link{PEdom.num}},
#' and \code{\link{PEdom.num1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2)
#' Tr<-rbind(A,B,C)
#' n<-10  #try also n<-20
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1,1,1)
#'
#' r<-1.4
#'
#' PEdom.num.tri(Xp,Tr,r,M)
#' IM<-inci.matPEtri(Xp,Tr,r,M)
#' dom.num.greedy #try also dom.num.exact(IM)
#'
#' gr.gam<-dom.num.greedy(IM)
#' gr.gam
#' Xp[gr.gam$i,]
#'
#' PEdom.num.tri(Xp,Tr,r,M=c(.4,.4))
#' }
#'
#' @export PEdom.num.tri
PEdom.num.tri <- function(Xp,tri,r,M=c(1,1,1))
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  n<-nrow(Xp)  #number of Xp points
  ind.tri<-mds<-mds.ind<-c()
  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],tri,boundary = TRUE)$i)
      ind.tri<-c(ind.tri,i) #indices of data points in the triangle, tri
  }

  Xtri<-matrix(Xp[ind.tri,],ncol=2) #data points in the triangle, tri

  ntri<-nrow(Xtri)  #number of points inside the triangle
  if (ntri==0)
  {gam<-0;
  res<-list(dom.num=gam, #domination number
            mds=NULL,  #a minimum dominating set
            ind.mds=NULL #indices of the mds
  )
  return(res); stop}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  Cl2e0<-cl2edgesMvert.reg(Xtri,tri,M)
  Cl2e<-Cl2e0$ext
  #for general r, points closest to opposite edges in the vertex regions
  Cl2e.ind<-Cl2e0$ind # indices of these extrema wrt Xtri
  Ext.ind =ind.tri[Cl2e.ind]
  #indices of these extrema wrt to the original data

  #Gamma=1 piece
  cnt<-0; j<-1;
  while (j<=3 & cnt==0)
  {
    if (Idom.num1PEtri(Cl2e[j,],Xtri,tri,r,M,rv=j)==1)
    {gam<-1; cnt<-1; mds<-rbind(mds,Cl2e[j,]);
    mds.ind=c(mds.ind,Ext.ind[j])
    } else
    {j<-j+1}
  }

  #Gamma=2 piece
  if (cnt==0)
  { k<-1; cnt2<-0;
  while (k<=2 & cnt2==0)
  {l<-k+1;
  while (l<=3 & cnt2==0)
  {
    if (Idom.num2PEtri(Cl2e[k,],Cl2e[l,],Xtri,tri,r,M,rv1=k,rv2=l)==1)
    {gam<-2;cnt2<-1 ; mds<-rbind(mds,Cl2e[c(k,l),]);
    mds.ind=c(mds.ind,Ext.ind[c(k,l)])
    } else {l<-l+1};
  }
  k<-k+1;
  }
  }

  #Gamma=3 piece
  if (cnt==0 && cnt2==0)
  {gam <-3; mds<-rbind(mds,Cl2e); mds.ind=c(mds.ind,Ext.ind)}

  row.names(mds)<-c()
  if (nrow(mds)==1){mds=as.vector(mds)}
  list(dom.num=gam, #domination number
       mds=mds, #a minimum dominating set
       ind.mds =mds.ind
       #indices of a minimum dominating set (wrt to original data)
  )
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point
#' in set \code{S} to the point \code{p} or
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p} in \eqn{N_{PE}(x,r)}
#' for some \eqn{x} in \code{S}\eqn{)}
#' for \code{S}, in the standard equilateral triangle,
#' that is, returns 1 if \code{p} is in \eqn{\cup_{x in S}N_{PE}(x,r)},
#' and returns 0 otherwise.
#'
#' PE proximity region is constructed
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' with the expansion parameter \eqn{r \ge 1}
#' and vertex regions are based
#' on center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)},
#' i.e., the center of mass of \eqn{T_e}
#' (which is equivalent to the circumcenter for \eqn{T_e}).
#'
#' Vertices of \eqn{T_e} are also labeled as 1, 2, and 3,
#' respectively.
#' If \code{p} is not in \code{S} and either \code{p}
#' or all points in \code{S} are outside \eqn{T_e}, it returns 0,
#' but if \code{p} is in \code{S},
#' then it always returns 1 regardless of its location
#' (i.e., loops are allowed).
#'
#' @param S A set of 2D points.
#' Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param p A 2D point.
#' Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region in the
#' standard equilateral triangle \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))};
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)}
#' i.e., the center of mass of \eqn{T_e}.
#'
#' @return \eqn{I(}\code{p} is in U_{x in \code{S}} \eqn{N_{PE}(x,r))}
#' for \code{S} in the standard equilateral triangle,
#' that is, returns 1 if \code{p} is in \code{S}
#' or inside \eqn{N_{PE}(x,r)} for at least
#' one \eqn{x} in \code{S}, and returns 0 otherwise.
#' PE proximity region is constructed with respect to the standard
#' equilateral triangle \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' with \code{M}-vertex regions
#'
#' @seealso \code{\link{IarcPEset2pnt.tri}}, \code{\link{IarcPEstd.tri}},
#' \code{\link{IarcPEtri}}, and \code{\link{IarcCSset2pnt.std.tri}}
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
#' r<-1.5
#'
#' S<-rbind(Xp[1,],Xp[2,])  #try also S<-c(.5,.5)
#' IarcPEset2pnt.std.tri(S,Xp[3,],r,M)
#' IarcPEset2pnt.std.tri(S,Xp[3,],r=1,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' IarcPEset2pnt.std.tri(S,Xp[3,],r,M)
#'
#' IarcPEset2pnt.std.tri(S,Xp[6,],r,M)
#' IarcPEset2pnt.std.tri(S,Xp[6,],r=1.25,M)
#'
#' P<-c(.4,.2)
#' S<-Xp[c(1,3,4),]
#' IarcPEset2pnt.std.tri(Xp,P,r,M)
#' }
#'
#' @export IarcPEset2pnt.std.tri
IarcPEset2pnt.std.tri <- function(S,p,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  if (dimension(M)==3)
  {M<-bary2cart(M,Te)}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (in.triangle(M,Te,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not in the interior of the triangle')}

  k<-nrow(S);
  dom<-0; i<-1;
  while (dom ==0 && i<= k)
  {
    if (IarcPEstd.tri(S[i,],p,r,M)==1)
    {dom<-1};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an arc from a point in set \code{S}
#' to the point \code{p} for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - one triangle case
#'
#' @description Returns \eqn{I(}\code{p} in \eqn{N_{PE}(x,r)}
#' for some \eqn{x} in \code{S}\eqn{)},
#' that is, returns 1 if \code{p} is in \eqn{\cup_{x in S}N_{PE}(x,r)},
#' and returns 0 otherwise.
#'
#' PE proximity region is constructed
#' with respect to the triangle \code{tri} with
#' the expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on the center, \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#'  default is \eqn{M=(1,1,1)}, i.e.,
#'  the center of mass of \code{tri}.
#' Vertices of \code{tri} are also labeled as 1, 2, and 3,
#' respectively.
#'
#' If \code{p} is not in \code{S} and either \code{p}
#' or all points in \code{S} are outside \code{tri}, it returns 0,
#' but if \code{p} is in \code{S},
#' then it always returns 1 regardless of its location
#' (i.e., loops are allowed).
#'
#' @param S A set of 2D points.
#' Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param p A 2D point.
#' Presence of an arc from a point in \code{S} to point \code{p} is checked
#' by the function.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region
#' constructed in the triangle \code{tri}; must be \eqn{\ge 1}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#'
#' @return I(\code{p} is in U_{x in \code{S}} N_{PE}(x,r)),
#' that is, returns 1 if \code{p} is in \code{S}
#' or inside \eqn{N_{PE}(x,r)} for at least
#' one \eqn{x} in \code{S}, and returns 0 otherwise,
#' where PE proximity region is constructed
#' with respect to the triangle \code{tri}
#'
#' @seealso \code{\link{IarcPEset2pnt.std.tri}}, \code{\link{IarcPEtri}},
#' \code{\link{IarcPEstd.tri}}, \code{\link{IarcASset2pnt.tri}},
#' and \code{\link{IarcCSset2pnt.tri}}
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
#' r<-1.5
#'
#' S<-rbind(Xp[1,],Xp[2,])  #try also S<-c(1.5,1)
#'
#' IarcPEset2pnt.tri(S,Xp[3,],Tr,r,M)
#' IarcPEset2pnt.tri(S,Xp[3,],r=1,Tr,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' IarcPEset2pnt.tri(S,Xp[3,],Tr,r,M)
#'
#' S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
#' IarcPEset2pnt.tri(S,Xp[3,],Tr,r,M)
#'
#' P<-c(.4,.2)
#' S<-Xp[c(1,3,4),]
#' IarcPEset2pnt.tri(Xp,P,Tr,r,M)
#' }
#'
#' @export IarcPEset2pnt.tri
IarcPEset2pnt.tri <- function(S,p,tri,r,M=c(1,1,1))
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  k<-nrow(S);
  dom<-0; i<-1;
  while (dom ==0 && i<= k)
  {
    if (IarcPEtri(S[i,],p,tri,r,M)==1)
    {dom<-1};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the set of points \code{S} being a dominating set
#' or not for Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{S} a dominating set of PE-PCD
#' whose vertices are the data points \code{Xp}\eqn{)}
#' for \code{S} in the standard equilateral triangle,
#' that is,
#' returns 1 if \code{S} is a dominating set of PE-PCD,
#' and returns 0 otherwise.
#'
#' PE proximity region is constructed
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(A,B,C)=T((0,0),(1,0),(1/2,\sqrt{3}/2))} with
#' expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_e}
#' (which is also equivalent to the circumcenter of \eqn{T_e}).
#' Vertices of \eqn{T_e} are also labeled as 1, 2, and 3,
#' respectively.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points
#' whose PE proximity regions are considered.
#' @param Xp A set of 2D points
#' which constitutes the vertices of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region in the
#' standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))};
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#'
#' @return \eqn{I(}\code{S} a dominating set of PE-PCD\eqn{)} for \code{S}
#' in the standard equilateral triangle,
#' that is, returns 1 if \code{S} is a dominating set of PE-PCD,
#' and returns 0 otherwise,
#' where PE proximity region is constructed in the standard equilateral triangle \eqn{T_e}.
#'
#' @seealso \code{\link{Idom.setPEtri}} and \code{\link{Idom.setCSstd.tri}}
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
#' r<-1.5
#'
#' S<-rbind(Xp[1,],Xp[2,])
#' Idom.setPEstd.tri(S,Xp,r,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,],c(.2,.5))
#' Idom.setPEstd.tri(S,Xp[3,],r,M)
#' }
#'
#' @export Idom.setPEstd.tri
Idom.setPEstd.tri <- function(S,Xp,r,M=c(1,1,1))
{
  if (!is.numeric(as.matrix(S)) ||
      !is.numeric(as.matrix(Xp)))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

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
    if (IarcPEset2pnt.std.tri(S,Xp[i,],r,M)==0)
      #this is where std equilateral triangle Te is implicitly used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The indicator for the set of points \code{S} being a dominating set
#' or not for Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description Returns \eqn{I(}\code{S} a dominating set of PE-PCD
#' whose vertices are the data set \code{Xp}\eqn{)}, that is,
#' returns 1 if \code{S} is a dominating set of PE-PCD,
#' and returns 0 otherwise.
#'
#' PE proximity region is constructed with
#' respect to the triangle \code{tri}
#' with the expansion parameter \eqn{r \ge 1} and vertex regions are based
#' on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' The triangle \code{tri}\eqn{=T(A,B,C)} has edges \eqn{AB}, \eqn{BC}, \eqn{AC}
#' which are also labeled as edges 3, 1, and 2, respectively.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param S A set of 2D points which is to be tested for
#' being a dominating set for the PE-PCDs.
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region
#' constructed in the triangle \code{tri}; must be \eqn{\ge 1}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#'
#' @return \eqn{I(}\code{S} a dominating set of PE-PCD\eqn{)},
#' that is, returns 1 if \code{S} is a dominating set of PE-PCD whose
#' vertices are the data points in \code{Xp}; and returns 0 otherwise,
#' where PE proximity region is constructed in
#' the triangle \code{tri}.
#'
#' @seealso \code{\link{Idom.setPEstd.tri}}, \code{\link{IarcPEset2pnt.tri}},
#' \code{\link{Idom.setCStri}}, and \code{\link{Idom.setAStri}}
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
#' r<-1.5
#'
#' S<-rbind(Xp[1,],Xp[2,])
#' Idom.setPEtri(S,Xp,Tr,r,M)
#'
#' S<-rbind(Xp[1,],Xp[2,],Xp[3,],Xp[5,])
#' Idom.setPEtri(S,Xp,Tr,r,M)
#'
#' S<-rbind(c(.1,.1),c(.3,.4),c(.5,.3))
#' Idom.setPEtri(S,Xp,Tr,r,M)
#' }
#'
#' @export Idom.setPEtri
Idom.setPEtri <- function(S,Xp,tri,r,M=c(1,1,1))
{
  if (!is.numeric(as.matrix(S)) ||
      !is.numeric(as.matrix(Xp)))
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

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  k<-nrow(S);
  n<-nrow(Xp);

  dom<-1; i<-1;
  while (dom ==1 && i<= n)
  {
    if (IarcPEset2pnt.tri(S,Xp[i,],tri,r,M)==0)
      #this is where tri is used
    {dom<-0};
    i<-i+1;
  }
  dom
} #end of the function
#'

#################################################################

#' @title The arcs of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data - one triangle case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends)
#' for data set \code{Xp} as the vertices of PE-PCD
#' and related parameters and the quantities of the digraph.
#'
#' PE proximity regions are constructed
#' with respect to the triangle \code{tri} with expansion
#' parameter \eqn{r \ge 1}, i.e.,
#' arcs may exist only for points inside \code{tri}.
#' It also provides various descriptions
#' and quantities about the arcs of the PE-PCD
#' such as number of arcs, arc density, etc.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of
#' the triangle \code{tri} or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#' \code{M}-vertex regions are recommended spatial inference,
#' due to geometry invariance property of the arc density
#' and domination number the PE-PCDs based on uniform data.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph,
#' the center \code{M} used to
#' construct the vertex regions and the expansion parameter \code{r}.}
#' \item{tess.points}{Points on which the tessellation of the study region
#' is performed, here, tessellation
#' is the support triangle.}
#' \item{tess.name}{Name of data set
#' (i.e. points from the non-target class) used in the tessellation
#' of the space (here, vertices of the triangle)}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitutes the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of PE-PCD
#' for 2D data set \code{Xp}
#' as vertices of the digraph}
#' \item{E}{Heads (or arrow ends) of the arcs of PE-PCD
#' for 2D data set \code{Xp}
#' as vertices of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph:
#' number of vertices, number of partition points,
#' number of triangles, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsPE}}, \code{\link{arcsAStri}},
#' and \code{\link{arcsCStri}}
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
#' r<-1.5  #try also r<-2
#'
#' Arcs<-arcsPEtri(Xp,Tr,r,M)
#' #or try with the default center Arcs<-arcsPEtri(Xp,Tr,r); M= (Arcs$param)$cent
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#'
#' #can add vertex regions
#' #but we first need to determine center is the circumcenter or not,
#' #see the description for more detail.
#' CC<-circumcenter.tri(Tr)
#' if (isTRUE(all.equal(M,CC)))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-prj.cent2edges(Tr,M)
#' }
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' #now we can add the vertex names and annotation
#' txt<-rbind(Tr,cent,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.03,-.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.03,.06,.04,.05,-.07)
#' txt.str<-c("A","B","C","M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export arcsPEtri
arcsPEtri <- function(Xp,tri,r,M=c(1,1,1))
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  n<-nrow(Xp)
  in.tri<-rep(0,n)
  for (i in 1:n)
    in.tri[i]<-in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri
  #indices the Xp points inside the triangle

  Xtri<-Xp[in.tri==1,] #the Xp points inside the triangle
  n2<-length(Xtri)/2

  #the arcs of PE-PCDs
  S<-E<-NULL #S is for source and E is for end points for the arcs
  if (n2>1)
  {
    for (j in 1:n2)
    {
      p1<-Xtri[j,];
      RV1<-ifelse(isTRUE(all.equal(M,CC)),
                  rel.vert.triCC(p1,tri)$rv,
                  rel.vert.tri(p1,tri,M)$rv) #vertex region for p1
      for (k in (1:n2)[-j])  #to avoid loops
      {
        p2<-Xtri[k,];
        if (IarcPEtri(p1,p2,tri,r,M,RV1)==1)
        {
          S <-rbind(S,p1); E <-rbind(E,p2);
        }
      }
    }
  }

  param<-list(M,r)
  Mr<-round(M,2)

  if (identical(M,"CC") || isTRUE(all.equal(M,CC)))
  {
    cname <-"CC"
    names(param)<-c("circumcenter","expansion parameter")
  main.txt<-paste("Arcs of PE-PCD \n with r = ",r," and Circumcenter",sep="")
  typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D Points in the Triangle with Expansion Parameter r = ",r," and Circumcenter",sep="")
  } else
  {
    cname <-"M"
    names(param)<-c("center","expansion parameter")
  main.txt<-paste("Arcs of PE-PCD\n with r = ",r," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
  typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D Points in the Triangle with Expansion Parameter r = ",r," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
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

#' @title The plot of the arcs of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for a 2D data set - one triangle case
#'
#' @description Plots the arcs of PE-PCD whose vertices are the data points, \code{Xp}
#' and the triangle \code{tri}.
#' PE proximity regions
#' are constructed with respect to the triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1}, i.e., arcs may exist only
#' for \code{Xp} points inside the triangle \code{tri}.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#' \code{M}-vertex regions are recommended spatial inference,
#' due to geometry invariance property of the arc density
#' and domination number the PE-PCDs based on uniform data.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param asp A \code{numeric} value,
#' giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by
#' typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes,
#' respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2,
#' giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param vert.reg A logical argument to add vertex regions to the plot,
#' default is \code{vert.reg=FALSE}.
#' @param \dots	Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of the PE-PCD
#' whose vertices are the points in data set \code{Xp}
#' and the triangle \code{tri}
#'
#' @seealso \code{\link{plotASarcs.tri}}, \code{\link{plotCSarcs.tri}},
#' and \code{\link{plotPEarcs}}
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
#' M<-as.numeric(runif.tri(1,Tr)$g)
#' #try also M<-c(1.6,1.0) or M<-circumcenter.tri(Tr)
#' r<-1.5  #try also r<-2
#' plotPEarcs.tri(Xp,Tr,r,M,main="Arcs of PE-PCD with r = 1.5",
#' xlab="",ylab="",vert.reg = TRUE)
#'
#' # or try the default center
#' #plotPEarcs.tri(Xp,Tr,r,main="Arcs of PE-PCD with r = 1.5",
#' #xlab="",ylab="",vert.reg = TRUE);
#' #M=(arcsPEtri(Xp,Tr,r)$param)$cent
#' #the part "M=(arcsPEtri(Xp,Tr,r)$param)$cent" is optional,
#' #for the below annotation of the plot
#'
#' #can add vertex labels and text to the figure (with vertex regions)
#' ifelse(isTRUE(all.equal(M,circumcenter.tri(Tr))),
#' {Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2); cent.name="CC"},
#' {Ds<-prj.cent2edges(Tr,M); cent.name="M"})
#'
#' txt<-rbind(Tr,M,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.04,-0.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.02,.07,.02,.04,-.06)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export plotPEarcs.tri
plotPEarcs.tri <- function(Xp,tri,r,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,
                           xlim=NULL,ylim=NULL,vert.reg=FALSE,...)
{
  arcsPE<-arcsPEtri(Xp,tri,r,M)
  S<-arcsPE$S
  E<-arcsPE$E
  cent = (arcsPE$param)$c

  Xp<-matrix(Xp,ncol=2)
  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1],cent[1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2],cent[2])}

  if ( isTRUE(all.equal(cent,circumcenter.tri(tri))) )
  {M="CC"}

  if (is.null(main))
  {if (identical(M,"CC")){
    main=paste("Arcs of PE-PCD\n with r = ",r," and Circumcenter",sep="")
  } else {Mr=round(cent,2)
  Mvec= paste(Mr, collapse=",")
  main=paste("Arcs of PE-PCD\n with r = ",r," and M = (",Mvec,")",sep="")}
  }

  if (vert.reg)
  {main=c(main,"\n (vertex regions added)")}

  plot(Xp,main=main,asp=asp, xlab=xlab, ylab=ylab,
       xlim=xlim,ylim=ylim,pch=".",cex=3,...)
  polygon(tri)
  if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}

  if (vert.reg){
    ifelse(isTRUE(all.equal(cent,circumcenter.tri(tri))),
           {A=tri[1,];B=tri[2,];C=tri[3,];
           Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2)},
           Ds<-prj.cent2edges(tri,M))
    L<-rbind(cent,cent,cent); R<-Ds
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }
} #end of the function
#'

#################################################################

#' @title The plot of the Proportional Edge (PE) Proximity Regions
#' for a 2D data set - one triangle case
#'
#' @description Plots the points in and outside of the triangle \code{tri}
#' and also the PE proximity regions
#' for points in data set \code{Xp}.
#'
#' PE proximity regions are defined
#' with respect to the triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1},
#' so PE proximity regions are defined only for points inside the
#' triangle \code{tri}.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#' \code{M}-vertex regions are recommended spatial inference,
#' due to geometry invariance property of the arc density
#' and domination number the PE-PCDs based on uniform data.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' for which PE proximity regions are constructed.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param asp A \code{numeric} value,
#' giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes,
#' respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2,
#' giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param vert.reg A logical argument to add vertex regions to the plot,
#' default is \code{vert.reg=FALSE}.
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the PE proximity regions for points
#' inside the triangle \code{tri}
#' (and just the points outside \code{tri})
#'
#' @seealso \code{\link{plotPEregs}}, \code{\link{plotASregs.tri}},
#' and \code{\link{plotCSregs.tri}}
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
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)
#' #try also M<-c(1.6,1.0) or M = circumcenter.tri(Tr)
#' r<-1.5  #try also r<-2
#'
#' plotPEregs.tri(Xp0,Tr,r,M)
#' Xp = Xp0[1,]
#' plotPEregs.tri(Xp,Tr,r,M)
#'
#' plotPEregs.tri(Xp,Tr,r,M,
#' main="PE Proximity Regions with r = 1.5",
#' xlab="",ylab="",vert.reg = TRUE)
#'
#' # or try the default center
#' #plotPEregs.tri(Xp,Tr,r,main="PE Proximity Regions with r = 1.5",xlab="",ylab="",vert.reg = TRUE);
#' #M=(arcsPEtri(Xp,Tr,r)$param)$c
#' #the part "M=(arcsPEtri(Xp,Tr,r)$param)$cent" is optional,
#' #for the below annotation of the plot
#'
#' #can add vertex labels and text to the figure (with vertex regions)
#' ifelse(isTRUE(all.equal(M,circumcenter.tri(Tr))),
#'        {Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2); cent.name="CC"},
#'        {Ds<-prj.cent2edges(Tr,M); cent.name<-"M"})
#'
#' txt<-rbind(Tr,M,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.03,-0.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.02,.07,.02,.05,-.06)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export plotPEregs.tri
plotPEregs.tri <- function(Xp,tri,r,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,vert.reg=FALSE,...)
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
  {stop('The triangle is degenerate')}

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!(is.point(M) || is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = circumcenter.tri(tri)
  if (isTRUE(all.equal(M,CC)))
  {cent=M; M="CC"}

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(identical(M,"CC") ||
        in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or
        not in the interior of the triangle')}

  n<-nrow(Xp)
  in.tri<-rep(0,n)
  for (i in 1:n)
    in.tri[i]<-in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri
  #indices of the Xp points inside the triangle

  Xtri<-matrix(Xp[in.tri==1,],ncol=2)
  #the Xp points inside the triangle
  nt<-length(Xtri)/2 #number of Xp points inside the triangle

ifelse(identical(M,"CC"),cent<-CC,cent <- M)
  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1],cent[1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2],cent[2])}

  xr<-xlim[2]-xlim[1]
  yr<-ylim[2]-ylim[1]

  Mr=round(cent,2)
  if (is.null(main))
  {if (identical(M,"CC")){
    main=paste("PE Proximity Regions\n with r = ",r," and Circumcenter",sep="")
  } else {Mvec= paste(Mr, collapse=",")
  main=paste("PE Proximity Regions\n with r = ",r," and M = (",Mvec,")",sep="")}
  }

  if (vert.reg)
  {main=c(main,"\n (vertex regions added)")}

  plot(Xp,main=main, asp=asp, xlab=xlab, ylab=ylab,xlim=xlim+xr*c(-.05,.05),
       ylim=ylim+yr*c(-.05,.05),pch=".",cex=3,...)
  polygon(tri,lty=2)
  if (nt>=1)
  {
    for (i in 1:nt)
    {
      P1<-Xtri[i,]
      RV<-ifelse(identical(M,"CC"),
                 rel.vert.triCC(P1,tri)$rv,
                 rel.vert.tri(P1,tri,M)$rv)

      pr<-NPEtri(P1,tri,r,M,rv=RV)
      polygon(pr,border="blue")
    }
  }

  if (vert.reg){
    ifelse(identical(M,"CC"),
           {A=tri[1,];B=tri[2,];C=tri[3,];
           Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2)},
           Ds<-prj.cent2edges(tri,cent))
    L<-rbind(cent,cent,cent); R<-Ds
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }

} #end of the function
#'

#################################################################

#' @title The arcs of Proportional Edge Proximity Catch Digraph (PE-PCD)
#' for 2D data - multiple triangle case
#'
#' @description
#' An object of class \code{"PCDs"}.
#' Returns arcs as tails (or sources) and heads (or arrow ends) of
#' Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data points in \code{Xp}
#' in the multiple triangle case
#' and related parameters and the quantities of the digraph.
#'
#' PE proximity regions are
#' defined with respect to the Delaunay triangles
#' based on \code{Yp} points with expansion parameter \eqn{r \ge 1} and
#' vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle or
#' based on circumcenter of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
#' Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' For the number of arcs, loops are not allowed so arcs are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the type of the digraph}
#' \item{parameters}{Parameters of the digraph,
#' the center used to construct the vertex regions
#' and the expansion parameter.}
#' \item{tess.points}{Points on which the tessellation
#' of the study region is performed,
#' here, tessellation
#' is Delaunay triangulation based on \code{Yp} points.}
#' \item{tess.name}{Name of data set used in tessellation,
#' it is \code{Yp} for this function}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitute the vertices of the digraph}
#' \item{S}{Tails (or sources) of the arcs of PE-PCD for 2D data set \code{Xp}
#' as vertices of the digraph}
#' \item{E}{Heads (or arrow ends) of the arcs of PE-PCD for 2D data set \code{Xp}
#' as vertices of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the digraph}
#' \item{quant}{Various quantities for the digraph: number of vertices,
#' number of partition points,
#' number of triangles, number of arcs, and arc density.}
#'
#' @seealso \code{\link{arcsPEtri}}, \code{\link{arcsAS}},
#' and \code{\link{arcsCS}}
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
#' r<-1.5  #try also r<-2
#'
#' Arcs<-arcsPE(Xp,Yp,r,M)
#' #or try with the default center Arcs<-arcsPE(Xp,Yp,r)
#' Arcs
#' summary(Arcs)
#' plot(Arcs)
#' }
#'
#' @export arcsPE
arcsPE <- function(Xp,Yp,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (nrow(Yp)==3)
  {
    res<-arcsPEtri(Xp,Yp,r,M)
  } else
  {
    if ((!is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates or
          "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)
    #the Xp points inside the convex hull of Yp

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
          Yi.tri<-as.basic.tri(Yi.Tri)$tri
          #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.basic.tri help page
          nl<-nrow(Xl)
          ifelse(identical(M,"CC"),
                 {rel.vert.ind<-rel.verts.triCC(Xl,Yi.tri)$rv;
          cent<-circumcenter.tri(Yi.tri)},
                 {rel.vert.ind<-rel.verts.tri(Xl,Yi.tri,M)$rv;
                 cent<-M})

          for (j in 1:nl)
          {RV<-rel.vert.ind[j]
          for (k in (1:nl)[-j])  # to avoid loops
            if (IarcPEtri(Xl[j,],Xl[k,],Yi.tri,r,cent,rv=RV)==1 )
            {
              S <-rbind(S,Xl[j,]); E <-rbind(E,Xl[k,]);
            }
          }
        }
      }
    }

    cname <-ifelse(identical(M,"CC"),"CC","M")
    param<-list(M,r)
    names(param)<-c("center","expansion parameter")

    if (identical(M,"CC")){
      main.txt=paste("Arcs of PE-PCD\n with r = ",r," and Circumcenter",sep="")
      typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D points in Multiple Triangles with Expansion Parameter r = ",r," and Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main.txt=paste("Arcs of PE-PCD\n with r = ",r," and Center ", cname," = (",Mvec,")",sep="")
    typ<-paste("Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D points in Multiple Triangles with Expansion parameter r = ",r," and Center ", cname," = (",Mvec,")",sep="")}

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

#' @title Incidence matrix for Proportional Edge Proximity Catch Digraphs
#' (PE-PCDs) - multiple triangle case
#'
#' @description Returns the incidence matrix of
#' Proportional Edge Proximity Catch Digraph (PE-PCD)
#' whose vertices are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' PE proximity regions are
#' defined with respect to the Delaunay triangles
#' based on \code{Yp} points with expansion parameter \eqn{r \ge 1} and
#' vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle
#' or based on circumcenter of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
#'
#' Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' For the incidence matrix loops are allowed,
#' so the diagonal entries are all equal to 1.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#'
#' @return Incidence matrix for the PE-PCD
#' with vertices being 2D data set, \code{Xp}.
#' PE proximity regions are constructed
#' with respect to the Delaunay triangles and \code{M}-vertex regions.
#'
#' @seealso \code{\link{inci.matPEtri}}, \code{\link{inci.matPEstd.tri}},
#' \code{\link{inci.matAS}}, and \code{\link{inci.matCS}}
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
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' r<-1.5  #try also r<-2
#'
#' IM<-inci.matPE(Xp,Yp,r,M)
#' IM
#' dom.num.greedy(IM)
#' #try also dom.num.exact(IM)
#' #might take a long time in this brute-force fashion ignoring the
#' #disconnected nature of the digraph inherent by the geometric construction of it
#' }
#'
#' @export inci.matPE
inci.matPE <- function(Xp,Yp,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (nrow(Yp)==3)
  {
    inc.mat<-inci.matPEtri(Xp,Yp,r,M)
  } else
  {
    if ((!is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates
          or "CC" for circumcenter')}

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
        Yi.tri<-as.basic.tri(Yi.Tri)$tri
        #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.basic.tri help page
        ifelse(identical(M,"CC"),{vert<-rel.vert.triCC(p1,Yi.tri)$rv;
        cent<-circumcenter.tri(Yi.tri)},
               {vert<-rel.vert.tri(p1,Yi.tri,M)$rv; cent<-M})

        for (j in 1:nx )
        {p2<-Xp[j,]
        inc.mat[i,j]<-IarcPEtri(p1,p2,Yi.tri,r,cent,rv=vert)
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

#' @title The plot of the arcs of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for a 2D data set - multiple triangle case
#'
#' @description Plots the arcs of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data
#' points in \code{Xp} in the multiple triangle case
#' and the Delaunay triangles based on \code{Yp} points.
#'
#' PE proximity regions are defined
#' with respect to the Delaunay triangles based on \code{Yp} points
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each Delaunay triangle
#' or based on circumcenter of
#' each Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
#' Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary
#' when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned by
#' the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' Loops are not allowed so arcs are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param asp A \code{numeric} value,
#' giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes,
#' respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2,
#' giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the arcs of the PE-PCD
#' whose vertices are the points in data set \code{Xp} and the Delaunay
#' triangles based on \code{Yp} points
#'
#' @seealso \code{\link{plotPEarcs.tri}}, \code{\link{plotASarcs}},
#' and \code{\link{plotCSarcs}}
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
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#'
#' r<-1.5  #try also r<-2
#'
#' plotPEarcs(Xp,Yp,r,M,xlab="",ylab="")
#' }
#'
#' @export plotPEarcs
plotPEarcs <- function(Xp,Yp,r,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
{
  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (nrow(Yp)==3)
  {
    plotPEarcs.tri(Xp,Yp,r,M,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    arcsPE<-arcsPE(Xp,Yp,r,M)
    S<-arcsPE$S
    E<-arcsPE$E

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    Xch<-Xin.convex.hullY(Xp,Yp)

    if (is.null(main))
    {if (identical(M,"CC")){
      main=paste("Arcs of PE-PCD\n with r = ",r," and Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main=paste("Arcs of PE-PCD\n with r = ",r," and M = (",Mvec,")",sep="")}
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
    plot(rbind(Xp),asp=asp,main=main, xlab=xlab, ylab=ylab,xlim=xlim,
         ylim=ylim,pch=".",cex=3,...)
    interp::plot.triSht(DTmesh, add=TRUE, do.points = TRUE)
    if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}
  }
} #end of the function
#'

#################################################################

#' @title The plot of the Proportional Edge (PE) Proximity Regions
#' for a 2D data set - multiple triangle case
#'
#' @description Plots the points in and outside of the Delaunay triangles
#' based on \code{Yp} points which partition
#' the convex hull of \code{Yp} points and also plots the PE proximity regions
#' for \code{Xp} points and the Delaunay triangles based on \code{Yp} points.
#'
#' PE proximity regions are constructed
#' with respect to the Delaunay triangles with the expansion parameter
#' \eqn{r \ge 1}.
#'
#' Vertex regions in each triangle is
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each Delaunay triangle
#' or based on circumcenter of
#' each Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:dom-num-NPE-Spat2011;textual}{pcds})
#' for more on the PE proximity regions.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' for which PE proximity regions are constructed.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#'  which is the center of mass of each triangle.
#' @param asp A \code{numeric} value,
#' giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes,
#' respectively (default=\code{NULL} for both)
#' @param xlim,ylim Two \code{numeric} vectors of length 2,
#' giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return Plot of the \code{Xp} points,
#' Delaunay triangles based on \code{Yp} points
#' and also the PE proximity regions
#' for \code{Xp} points inside the convex hull of \code{Yp} points
#'
#' @seealso \code{\link{plotPEregs.tri}}, \code{\link{plotASregs}},
#' and \code{\link{plotCSregs}}
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
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#' r<-1.5  #try also r<-2
#'
#' plotPEregs(Xp,Yp,r,M,xlab="",ylab="")
#' }
#'
#' @export plotPEregs
plotPEregs <- function(Xp,Yp,r,M=c(1,1,1),asp=NA,main=NULL,xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,...)
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (nrow(Yp)==3)
  {
    plotPEregs.tri(Xp,Yp,r,M,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    if ((!is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates
          or "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)
    #the Xp points inside the convex hull of Yp points

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)  #number of Delaunay triangles
    nx2<-nrow(Xch)  #number of Xp points inside the convex hull of Yp points

    if (nx2>=1)
    {
      i.tr<-rep(0,nx2)
      #the vector of indices for the triangles that contain the Xp points
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
    xlim<-xlim+xr*c(-.05,.05)
    }
    if (is.null(ylim))
    {ylim<-range(Yp[,2],Xp[,2])
    yr<-ylim[2]-ylim[1]
    ylim<-ylim+yr*c(-.05,.05)
    }

    if (is.null(main))
    {if (identical(M,"CC")) {
      main=paste("PE Proximity Regions\n with r = ",r," and Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main=paste("PE Proximity Regions\n with r = ",r," and M = (",Mvec,")",sep="")}
    }

    plot(rbind(Xp),asp=asp,main=main, xlab=xlab, ylab=ylab,
         xlim=xlim,ylim=ylim,pch=".",cex=3,...)

    for (i in 1:nt)
    {
      Tri<-Yp[DTr[i,],]  #vertices of the ith triangle
      tri<-as.basic.tri(Tri)$tri
      #convert the triangle Tri into an nonscaled basic triangle, see as.basic.tri help page

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
            ifelse(identical(M,"CC"),{RV<-rel.vert.triCC(P1,tri)$rv; cent<-circumcenter.tri(tri)},
                   {RV<-rel.vert.tri(P1,tri,M)$rv; cent<-M})

            pr<-NPEtri(P1,tri,r,cent,rv=RV)
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

#' @title The domination number of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) - multiple triangle case
#'
#' @description Returns the domination number,
#' indices of a minimum dominating set of PE-PCD whose vertices are the data
#' points in \code{Xp} in the multiple triangle case
#' and domination numbers for the Delaunay triangles
#' based on \code{Yp} points.
#'
#' PE proximity regions are defined
#' with respect to the Delaunay triangles based on \code{Yp} points
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle or based on
#' circumcenter of each Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the
#' triangle). Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the
#' same type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' Loops are allowed for the domination number.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds})
#' for more on the domination number of PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and
#' the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#'  (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#'
#' @return A \code{list} with three elements
#' \item{dom.num}{Domination number of the PE-PCD
#' whose vertices are \code{Xp} points.
#' PE proximity regions are constructed
#' with respect to the Delaunay triangles
#' based on the \code{Yp} points with expansion parameter \eqn{r \ge 1}.}
#' #\item{mds}{A minimum dominating set of the PE-PCD
#' whose vertices are \code{Xp} points}
#' \item{ind.mds}{The vector of data indices of the minimum dominating set
#' of the PE-PCD whose vertices are \code{Xp} points.}
#' \item{tri.dom.nums}{The vector of domination numbers
#' of the PE-PCD components
#' for the Delaunay triangles.}
#'
#' @seealso \code{\link{PEdom.num.tri}}, \code{\link{PEdom.num.tetra}},
#' \code{\link{dom.num.exact}}, and \code{\link{dom.num.greedy}}
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
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' M<-c(1,1,1)  #try also M<-c(1,2,3)
#' r<-1.5  #try also r<-2
#' PEdom.num(Xp,Yp,r,M)
#' }
#'
#' @export PEdom.num
PEdom.num <- function(Xp,Yp,r,M=c(1,1,1))
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if ((!is.point(M,3) && M!="CC") || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates
        or "CC" for circumcenter')}

  if (nrow(Yp)==3)
  {
    res<-PEdom.num.tri(Xp,Yp,r,M)
  } else
  {
    n<-nrow(Xp)  #number of Xp points
    m<-nrow(Yp)  #number of Yp points

    Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    #Delaunay triangulation
    Ytri<-matrix(interp::triangles(Ytrimesh)[,1:3],ncol=3);
    #the Delaunay triangles
    nt<-nrow(Ytri)  #number of Delaunay triangles
    inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
    #logical indices for Xp points in convex hull of Yp points
    Xch<-matrix(Xp[inCH==TRUE,],ncol=2)
    indCH=which(inCH) #vector of indices of the data in the convex hull

    gam<-rep(0,nt);  mds<-mds.ind<-c()
    if (nrow(Xch)>=1)
    {
      Tri.Ind<-indices.delaunay.tri(Xch,Yp,Ytrimesh)
      #indices of triangles in which the points in the data fall

      #calculation of the domination number
      for (i in 1:nt)
      {
        ith.tri.ind = Tri.Ind==i
      Xpi<-matrix(Xch[ith.tri.ind,],ncol=2)  #points in ith Delaunay triangle
      indCHi = indCH[ith.tri.ind]
      #indices of Xpi points (wrt to original data indices)

        ni<-nrow(Xpi)  #number of points in ith triangle
        if (ni==0)
        {
          gam[i]<-0
        } else
        {
          Yi.Tri<-Yp[Ytri[i,],] #vertices of ith triangle
          Yi.tri<-as.basic.tri(Yi.Tri)$tri
          ifelse(identical(M,"CC"),
                 {cent<-circumcenter.tri(Yi.tri);
                 cl2v <- cl2edgesCCvert.reg(Xpi,Yi.tri);
                 Clvert<-cl2v$ext;
                 Clvert.ind<-cl2v$ind # indices of these extrema wrt Xpi
                 },
                 {cent<-M;
                 cl2v <- cl2edgesMvert.reg(Xpi,Yi.tri,cent);
                 Clvert<-cl2v$ext;
                 Clvert.ind<-cl2v$ind # indices of these extrema wrt Xpi}
                 })

          Ext.ind = indCHi[Clvert.ind]
          #indices of these extrema wrt to the original data

          #Gamma=1 piece
          cnt<-0; j<-1;
          while (j<=3 & cnt==0)
          {
            if (Idom.num1PEtri(Clvert[j,],Xpi,Yi.tri,r,cent,rv=j)==1)
            {gam[i]<-1; cnt<-1;
            mds<-rbind(mds,Clvert[j,]);
            mds.ind=c(mds.ind,Ext.ind[j])
            } else
            {j<-j+1}
          }

          #Gamma=2 piece
          if (cnt==0)
          { k<-1; cnt2<-0;
          while (k<=2 & cnt2==0)
          {l<-k+1;
          while (l<=3 & cnt2==0)
          {
            if (Idom.num2PEtri(Clvert[k,],Clvert[l,],Xpi,Yi.tri,r,cent,
                          rv1=k,rv2=l)==1)
            {gam[i]<-2;cnt2<-1; mds<-rbind(mds,Clvert[c(k,l),]);
            mds.ind=c(mds.ind,Ext.ind[c(k,l)])
            } else {l<-l+1};
          }
          k<-k+1;
          }
          }

          if (cnt==0 && cnt2==0)
          {gam[i]<-3; mds<-rbind(mds,Clvert); mds.ind=c(mds.ind,Ext.ind)}
        }

      }
    }

    Gam<-sum(gam)  #domination number for the entire digraph
    row.names(mds)<-c()

    res<-list(dom.num=Gam, #domination number
            #  mds=mds, #a minimum dominating set
            ind.mds =mds.ind,
            #indices of a minimum dominating set (wrt to original data)
              tri.dom.nums=gam #domination numbers for the Delaunay triangles
    )
    }

res
} #end of the function
#'

#################################################################

#' @title The domination number of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) with non-degeneracy centers - multiple triangle case
#'
#' @description Returns the domination number,
#' indices of a minimum dominating set of PE-PCD
#' whose vertices are the data
#' points in \code{Xp} in the multiple triangle case
#' and domination numbers for the Delaunay triangles based on \code{Yp} points
#' when PE-PCD is constructed with vertex regions
#' based on non-degeneracy centers.
#'
#' PE proximity regions are defined
#' with respect to the Delaunay triangles based on \code{Yp} points
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions in each triangle are
#' based on the center \eqn{M}
#' which is one of the 3 centers
#' that renders the asymptotic distribution of domination number
#' to be non-degenerate for a given value of \code{r} in \eqn{(1,1.5)}
#' and \code{M} is center of mass for \eqn{r=1.5}.
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' Loops are allowed for the domination number.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011,ceyhan:mcap2012;textual}{pcds})
#' more on the domination number of PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds})
#' for more on Delaunay triangulation and
#' the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]} here.
#'
#' @return A \code{list} with three elements
#' \item{dom.num}{Domination number of the PE-PCD
#' whose vertices are \code{Xp} points. PE proximity regions are
#' constructed with respect to the Delaunay triangles
#' based on the \code{Yp} points with expansion parameter \eqn{r in (1,1.5]}.}
#' #\item{mds}{A minimum dominating set of the PE-PCD
#' whose vertices are \code{Xp} points.}
#' \item{ind.mds}{The data indices of the minimum dominating set of the PE-PCD
#' whose vertices are \code{Xp} points.}
#' \item{tri.dom.nums}{Domination numbers of the PE-PCD components
#' for the Delaunay triangles.}
#'
#' @seealso \code{\link{PEdom.num.tri}}, \code{\link{PEdom.num.tetra}},
#' \code{\link{dom.num.exact}}, and \code{\link{dom.num.greedy}}
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
#' r<-1.5  #try also r<-2
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' PEdom.num.nondeg(Xp,Yp,r)
#' }
#'
#' @export PEdom.num.nondeg
PEdom.num.nondeg <- function(Xp,Yp,r)
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

  if (!is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (nrow(Yp)==3)
  { rcent<-sample(1:3,1)  #random center selection from M1,M2,M3
  cent.nd<-center.nondegPE(Yp,r)[rcent,]
  res<-PEdom.num.tri(Xp,Yp,r,cent.nd)
  } else
  {
    Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    #Delaunay triangulation
    Ytri<-matrix(interp::triangles(Ytrimesh)[,1:3],ncol=3);
    #the Delaunay triangles
    nt<-nrow(Ytri)  #number of Delaunay triangles
    inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
    #logical indices for Xp points in convex hull of Yp points
    Xch<-matrix(Xp[inCH==TRUE,],ncol=2)
    indCH=which(inCH) #vector of indices of the data in the convex hull

    gam<-rep(0,nt);  mds<-mds.ind<-c()
    if (nrow(Xch)>=1)
    {
      Tri.Ind<-indices.delaunay.tri(Xch,Yp,Ytrimesh)
      #indices of triangles in which the points in the data fall

      #calculation of the domination number
      for (i in 1:nt)
      { ith.tri.ind = Tri.Ind==i
        Xpi<-matrix(Xch[ith.tri.ind,],ncol=2)
        #points in ith Delaunay triangle
        indCHi = indCH[ith.tri.ind]
        #indices of Xpi points (wrt to original data indices)

        ni<-nrow(Xpi)  #number of points in ith triangle
        if (ni==0)
        {
          gam[i]<-0
        } else
        {
          Yi.tri<-Yp[Ytri[i,],] #vertices of ith triangle
          rcent<-sample(1:3,1)  #random center selection from M1,M2,M3
          Centi<-center.nondegPE(Yi.tri,r)[rcent,]
          cl2v = cl2edgesMvert.reg(Xpi,Yi.tri,Centi)
          Clvert<-cl2v$ext
          #for general r, points closest to opposite edges in the vertex regions
          Clvert.ind<-cl2v$ind # indices of these extrema wrt Xpi
          Ext.ind =indCHi[Clvert.ind]
          #indices of these extrema wrt to the original data

          #Gamma=1 piece
          cnt<-0; j<-1;
          while (j<=3 & cnt==0)
          {
            if (Idom.num1PEtri(Clvert[j,],Xpi,Yi.tri,r,Centi,rv=j)==1)
            {gam[i]<-1; cnt<-1; mds<-rbind(mds,Clvert[j,]);
            mds.ind=c(mds.ind,Ext.ind[j])
            } else
            {j<-j+1}
          }

          #Gamma=2 piece
          if (cnt==0)
          { k<-1; cnt2<-0;
          while (k<=2 & cnt2==0)
          {l<-k+1;
          while (l<=3 & cnt2==0)
          {
            if (Idom.num2PEtri(Clvert[k,],Clvert[l,],Xpi,Yi.tri,r,Centi,
                          rv1=k,rv2=l)==1)
            {gam[i]<-2;cnt2<-1; mds<-rbind(mds,Clvert[c(k,l),]);
            mds.ind=c(mds.ind,Ext.ind[c(k,l)])
            } else {l<-l+1};
          }
          k<-k+1;
          }
          }

          if (cnt==0 && cnt2==0)
          {gam[i]<-3; mds<-rbind(mds,Clvert); mds.ind=c(mds.ind,Ext.ind)}

        }

      }
    }

    Gam<-sum(gam)  #domination number for the entire digraph
    row.names(mds)<-c()

    res<-list(dom.num=Gam, #domination number
             # mds=mds, #a minimum dominating set
              ind.mds =mds.ind,
             #indices of a minimum dominating set (wrt to original data)
              tri.dom.nums=gam #domination numbers for the Delaunay triangles
    )
  }
  res
} #end of the function
#'

#################################################################

#' @title Asymptotic probability that domination number of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) equals 2
#' where vertices of the digraph are uniform points in a triangle
#'
#' @description Returns \eqn{P(}domination number\eqn{=2)}
#' for PE-PCD for uniform data in a triangle,
#' when the sample size \eqn{n} goes to
#' infinity (i.e., asymptotic probability of domination number \eqn{= 2}).
#'
#' PE proximity regions are constructed
#' with respect to the triangle
#' with the expansion parameter \eqn{r \ge 1} and
#' \eqn{M}-vertex regions where \eqn{M} is the vertex
#' that renders the asymptotic distribution of the domination
#' number non-degenerate for the given value of \code{r} in \eqn{(1,1.5]}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]} to attain non-degenerate asymptotic distribution
#' for the domination number.
#'
#' @return \eqn{P(}domination number\eqn{=2)}
#' for PE-PCD for uniform data on an triangle as the sample size \eqn{n}
#' goes to infinity
#'
#' @seealso \code{\link{Pdom.num2PE1D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' Pdom.num2PEtri(r=1.5)
#' Pdom.num2PEtri(r=1.4999999999)
#'
#' Pdom.num2PEtri(r=1.5) / Pdom.num2PEtri(r=1.4999999999)
#'
#' rseq<-seq(1.01,1.49999999999,l=20)  #try also l=100
#' lrseq<-length(rseq)
#'
#' pg2<-vector()
#' for (i in 1:lrseq)
#' {
#'   pg2<-c(pg2,Pdom.num2PEtri(rseq[i]))
#' }
#'
#' plot(rseq, pg2,type="l",xlab="r",
#' ylab=expression(paste("P(", gamma, "=2)")),
#'      lty=1,xlim=range(rseq)+c(0,.01),ylim=c(0,1))
#' points(rbind(c(1.50,Pdom.num2PEtri(1.50))),pch=".",cex=3)
#' }
#'
#' @export Pdom.num2PEtri
Pdom.num2PEtri <- function(r)
{
  if (!is.point(r,1) || r<=1 || r>1.5)
  {stop('the argument must be a scalar in (1,1.5]')}

  if (r==1.5)
  {pg2<-.7413}
  else
  {pg2<--(1/2)*(pi*r^2-2*atan(r*(r-1)/sqrt(-r^4+2*r^3-r^2+1))*r^2-pi*r+2*atan(r*(r-1)/sqrt(-r^4+2*r^3-r^2+1))*r-2*sqrt(-r^4+2*r^3-r^2+1))/(-r^4+2*r^3-r^2+1)^(3/2)}
  pg2
} #end of the function
#'

#################################################################

#' @title A test of segregation/association based on domination number of
#' Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D data -
#' Binomial Approximation
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function
#' which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster away from \code{Yp} points
#' i.e., cluster around the centers of the Delaunay triangles)
#' and association (where \code{Xp} points cluster around \code{Yp} points)
#' based on the (asymptotic) binomial distribution of the
#' domination number of PE-PCD for uniform 2D data
#' in the convex hull of \code{Yp} points.
#'
#' The function yields the test statistic,
#' \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval,
#' estimate and null value for the parameter of interest
#' (which is \eqn{Pr(}domination number\eqn{\le 2)}),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points, probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{\le 2)}) equals
#' to its expected value under the uniform distribution) and
#' \code{alternative} could be two-sided, or right-sided
#' (i.e., data is accumulated
#' around the \code{Yp} points, or association)
#' or left-sided (i.e., data is accumulated
#' around the centers of the triangles, or segregation).
#'
#' PE proximity region is constructed
#' with the expansion parameter \eqn{r \ge 1} and \eqn{M}-vertex regions
#' where \eqn{M} is a center
#' that yields non-degenerate asymptotic distribution
#' of the domination number.
#'
#' The test statistic is based on the binomial distribution,
#' when success is defined as domination number being less than
#' or equal to 2 in the one triangle case
#' (i.e., number of failures is equal
#' to number of times restricted domination number = 3
#' in the triangles).
#' That is, the test statistic is based on the domination number
#' for \code{Xp} points inside convex hull of \code{Yp} points
#' for the PE-PCD and default convex hull correction, \code{ch.cor},
#' is \code{FALSE} where \code{M} is the center
#' that yields nondegenerate asymptotic distribution
#' for the domination number.
#' For this approximation to work,
#' number of \code{Xp} points must be at least 7 times more than
#' number of \code{Yp} points.
#'
#' PE proximity region is constructed
#' with the expansion parameter \eqn{r \ge 1} and \eqn{CM}-vertex regions
#' (i.e., the test is not available for a general center \eqn{M}
#' at this version of the function).
#'
#' **Caveat:** This test is currently a conditional test,
#' where \code{Xp} points are assumed to be random,
#' while \code{Yp} points are
#' assumed to be fixed (i.e., the test is conditional on \code{Yp} points).
#' Furthermore, the test is a large sample test
#' when \code{Xp} points are substantially larger than \code{Yp} points,
#' say at least 7 times more.
#' This test is more appropriate
#' when supports of \code{Xp} and \code{Yp} have a substantial overlap.
#' Currently, the \code{Xp} points
#' outside the convex hull of \code{Yp} points
#' are handled with a convex hull correction factor
#' (see the description below and the function code.)
#' Removing the conditioning
#' and extending it to the case of non-concurring supports is
#' an ongoing topic of research of the author of the package.
#'
#' See also (\insertCite{ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]}.
#' @param ch.cor A logical argument for convex hull correction,
#' default \code{ch.cor=FALSE},
#' recommended when both \code{Xp} and \code{Yp}
#' have the same rectangular support.
#' @param ndt Number of Delaunay triangles based on \code{Yp} points,
#' default is \code{NULL}.
#' @param alternative Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval,
#' default is \code{0.95}, for the probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{=3)} for PE-PCD
#' whose vertices are the 2D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test
#' for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval
#' for \eqn{Pr(}Domination Number\eqn{\le 2)}
#' at the given level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{A \code{vector} with two entries:
#' first is is the estimate of the parameter, i.e.,
#' \eqn{Pr(}Domination Number\eqn{=3)} and second is the domination number}
#' \item{null.value}{Hypothesized value for the parameter,
#' i.e., the null value for \eqn{Pr(}Domination Number\eqn{\le 2)}}
#' \item{alternative}{Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{PEdom.num.norm.test}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-100; ny<-5 #try also nx<-1000; ny<-10
#' r<-1.4  #try also r<-1.5
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' plotDelaunay.tri(Xp,Yp,xlab="",ylab="")
#' PEdom.num.binom.test(Xp,Yp,r) #try also #PEdom.num.binom.test(Xp,Yp,r,alt="l") and
#' # PEdom.num.binom.test(Xp,Yp,r,alt="g")
#' PEdom.num.binom.test(Xp,Yp,r,ch=TRUE)
#'
#' #or try
#' ndt<-num.delaunay.tri(Yp)
#' PEdom.num.binom.test(Xp,Yp,r,ndt=ndt)
#' #values might differ due to the random of choice of the three centers M1,M2,M3
#' #for the non-degenerate asymptotic distribution of the domination number
#' }
#'
#' @export PEdom.num.binom.test
PEdom.num.binom.test <- function(Xp,Yp,r,ch.cor=FALSE,ndt=NULL,alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one of \"greater\", \"less\", \"two.sided\"")

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

  if (!is.point(r,1) || r<=1 || r>1.5)
  {stop('r must be in (1,1.5] for domination number to be asymptotically non-degenerate')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) ||
        conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  if (is.null(ndt))
  {ndt<-num.delaunay.tri(Yp)} #number of Delaunay triangles

  Dom.num = PEdom.num.nondeg(Xp,Yp,r)
  Gammas<-Dom.num$tri.dom
  #vector of domination numbers for the Delaunay triangles

  estimate1<-Dom.num$dom.num #domination number of the entire PE-PCD
 # ind0<- Gammas>0; Gammas=(3-Gammas0)
  Bm<-sum(Gammas<=2)
  #sum((3-Gammas)[ind0]>0) #sum(Gammas-2>0);
  #the binomial test statistic, success is dom num <= 2

  method <-c("Large Sample Binomial Test based on the Domination Number of PE-PCD for Testing Uniformity of 2D Data ---")

  if (ch.cor==TRUE)  #the part for the convex hull correction
  {
    nx<-nrow(Xp)  #number of Xp points
    ny<-nrow(Yp)  #number of Yp points
    Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    #Delaunay triangulation

    inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
    #logical indices for Xp points in convex hull of Yp points

    outch<-nx-sum(inCH)
    prop.out<-outch/nx #observed proportion of points outside convex hull
    exp.prop.out<-1.7932/ny+1.2229/sqrt(ny)
    #expected proportion of points outside convex hull

    Bm<-Bm*(1-(prop.out-exp.prop.out))
    method <-c(method, " with Convex Hull Correction")
  } else
  {    method <-c(method, " without Convex Hull Correction")}

  p<-Pdom.num2PEtri(r)
  x<-round(Bm)
  pval <-switch(alternative, less = pbinom(x, ndt, p),
                greater = pbinom(x - 1, ndt, p, lower.tail = FALSE),
                two.sided = {if (p == 0) (x == 0) else if (p == 1) (x == ndt)
                  else { relErr <-1 + 1e-07
                  d <-dbinom(x, ndt, p)
                  m <-ndt * p
                  if (x == m) 1 else if (x < m)
                  {i <-seq.int(from = ceiling(m), to = ndt)
                  y <-sum(dbinom(i, ndt, p) <= d * relErr)
                  pbinom(x, ndt, p) + pbinom(ndt - y, ndt, p, lower.tail = FALSE)
                  } else {
                    i <-seq.int(from = 0, to = floor(m))
                    y <-sum(dbinom(i, ndt, p) <= d * relErr)
                    pbinom(y - 1, ndt, p) + pbinom(x - 1, ndt, p, lower.tail = FALSE)
                  }
                  }
                })

  p.L <- function(x, alpha) {
    if (x == 0)
      0
    else qbeta(alpha, x, ndt - x + 1)
  }
  p.U <- function(x, alpha) {
    if (x == ndt)
      1
    else qbeta(1 - alpha, x + 1, ndt - x)
  }
  cint <-switch(alternative, less = c(0, p.U(x, 1 - conf.level)),
                greater = c(p.L(x, 1 - conf.level), 1), two.sided = {
                  alpha <-(1 - conf.level)/2
                  c(p.L(x, alpha), p.U(x, alpha))
                })
  attr(cint, "conf.level") <- conf.level

  estimate2 <-x/ndt
  names(x) <- "#(domination number is <= 2)"
  #"domination number - 2 * number of Delaunay triangles"
  names(ndt) <-"number of Delaunay triangles based on Yp points"
  #paste("number of Delaunay triangles based on ", yname,sep="")

  names(p) <-"Pr(Domination Number <=2)"

  names(estimate1) <-c(" domination number  ")
  names(estimate2) <-c("|| Pr(domination number <= 2)")

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

#' @title A test of segregation/association based on domination number of
#' Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D data -
#' Normal Approximation
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function
#' which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster
#' away from \code{Yp} points i.e.,
#' cluster around the centers of the Delaunay
#' triangles) and association (where \code{Xp} points cluster
#' around \code{Yp} points) based on the normal approximation
#' to the binomial distribution of the domination number of PE-PCD
#' for uniform 2D data
#' in the convex hull of \code{Yp} points
#'
#' The function yields the test statistic, \eqn{p}-value
#' for the corresponding \code{alternative},
#' the confidence interval, estimate and null value
#' for the parameter of interest
#' (which is \eqn{Pr(}domination number\eqn{\le 2)}),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points, probability of success
#' (i.e., \eqn{Pr(}domination number\eqn{\le 2)}) equals
#' to its expected value under the uniform distribution) and
#' \code{alternative} could be two-sided, or right-sided
#' (i.e., data is accumulated around the \code{Yp} points, or association)
#' or left-sided (i.e., data is accumulated
#' around the centers of the triangles,
#' or segregation).
#'
#' PE proximity region is constructed
#' with the expansion parameter \eqn{r \ge 1}
#' and \eqn{M}-vertex regions where M
#' is a center that yields non-degenerate asymptotic distribution of
#' the domination number.
#'
#' The test statistic is based on the normal approximation
#' to the binomial distribution,
#' when success is defined as domination number being less than
#' or equal to 2 in the one triangle case
#' (i.e., number of failures is equal to number of times
#' restricted domination number = 3
#' in the triangles).
#' That is, the test statistic is
#' based on the domination number for \code{Xp} points
#' inside convex hull of \code{Yp}
#' points for the PE-PCD and default convex hull correction, \code{ch.cor},
#' is \code{FALSE}
#' where \code{M} is the center
#' that yields nondegenerate asymptotic distribution
#' for the domination number.
#'
#' For this approximation to work,
#' number of \code{Yp} points must be at least 5
#' (i.e., about 7 or more Delaunay triangles)
#' and number of \code{Xp} points must be at least 7 times more than
#' the number of \code{Yp} points.
#'
#' See also (\insertCite{ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be in \eqn{(1,1.5]}.
#' @param ch.cor A logical argument for convex hull correction,
#' default \code{ch.cor=FALSE},
#' recommended when both \code{Xp} and \code{Yp}
#' have the same rectangular support.
#' @param ndt Number of Delaunay triangles based on \code{Yp} points,
#' default is \code{NULL}.
#' @param alternative Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval,
#' default is \code{0.95}, for the domination number of
#' PE-PCD whose vertices are the 2D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test
#' for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval for the domination number
#' at the given level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{A \code{vector} with two entries:
#' first is the domination number,
#' and second is the estimate of the parameter, i.e.,
#' \eqn{Pr(}Domination Number\eqn{=3)}}
#' \item{null.value}{Hypothesized value for the parameter,
#' i.e., the null value for expected domination number}
#' \item{alternative}{Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{PEdom.num.binom.test}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' nx<-100; ny<-5 #try also nx<-1000; ny<-10
#' r<-1.5  #try also r<-2 or r<-1.25
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#' #try also Yp<-cbind(runif(ny,0,1),runif(ny,0,1))
#'
#' plotDelaunay.tri(Xp,Yp,xlab="",ylab="")
#' PEdom.num.norm.test(Xp,Yp,r) #try also PEdom.num.norm.test(Xp,Yp,r, alt="l")
#'
#' PEdom.num.norm.test(Xp,Yp,1.25,ch=TRUE)
#'
#' #or try
#' ndt<-num.delaunay.tri(Yp)
#' PEdom.num.norm.test(Xp,Yp,r,ndt=ndt)
#' #values might differ due to the random of choice of the three centers M1,M2,M3
#' #for the non-degenerate asymptotic distribution of the domination number
#' }
#'
#' @export PEdom.num.norm.test
PEdom.num.norm.test <- function(Xp,Yp,r,ch.cor=FALSE,ndt=NULL,alternative=c("two.sided", "less", "greater"),conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one of \"greater\", \"less\", \"two.sided\"")

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

  if (!is.point(r,1) || r<=1 || r>1.5)
  {stop('r must be in (1,1.5] for domination number to be asymptotically non-degenerate')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) ||
        conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  if (is.null(ndt))
  {ndt<-num.delaunay.tri(Yp)} #number of Delaunay triangles

  Dom.num = PEdom.num.nondeg(Xp,Yp,r)
  Gammas<-Dom.num$tri.dom #domination numbers for the Delaunay triangles

  estimate1<-Gam<-Dom.num$dom.num #domination number of the entire PE-PCD
 # ind0<- Gammas>0;  Gammas=(3-Gammas0)
  Bm<-sum(Gammas<=2) #sum((3-Gammas)[ind0]>0) #sum(Gammas-2>0);
  #the binomial test statistic, success is dom num <= 2

  p<-Pdom.num2PEtri(r)
  Exp.Gam <-ndt*p  #expected Bm value
  estimate2<-Bm/ndt #estimated Pr(gamma<=2)
  st.err<-sqrt(ndt*p*(1-p))
  TS0<-(Bm-Exp.Gam)/st.err  #the standardized test statistic
  #TS0<-(Gam-Exp.Gam)/st.err  #the standardized test statistic

  method <-c("Normal Approximation to the Domination Number of PE-PCD for Testing Uniformity of 2D Data ---")
  if (ch.cor==FALSE)  #the part for the convex hull correction
  {
    TS<-TS0
    method <-c(method, " without Convex Hull Correction")
  } else
  {
    n<-nrow(Xp)  #number of Xp points
    m<-nrow(Yp)  #number of Yp points
    Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    #Delaunay triangulation

    inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
    #logical indices for Xp points in convex hull of Yp points

    outch<-n-sum(inCH)
    prop.out<-outch/n #observed proportion of points outside convex hull
    exp.prop.out<-1.7932/m+1.2229/sqrt(m)
    #expected proportion of points outside convex hull

    TS<-TS0*(1-(prop.out-exp.prop.out))
    method <-c(method, " with Convex Hull Correction")
  }

  names(estimate1) <-c(" domination number  ")
  names(estimate2) <-c("|| Pr(domination number <= 2)")

  null.gam<-Exp.Gam
  names(null.gam) <-c("expected domination number")
  names(TS) <-"standardized domination number (i.e., Z)"

  if (alternative == "less") {
    pval <-pnorm(TS)
    cint <-Gam+c(-Inf, qnorm(conf.level))*st.err
  }
  else if (alternative == "greater") {
    pval <-pnorm(TS, lower.tail = FALSE)
    cint <-Gam+c(-qnorm(conf.level),Inf)*st.err
  }
  else {
    pval <-2 * pnorm(-abs(TS))
    alpha <-1 - conf.level
    cint <-qnorm(1 - alpha/2)
    cint <-Gam+c(-cint, cint)*st.err
  }

  attr(cint, "conf.level") <- conf.level

  rval <-list(
    statistic=TS,
    p.value=pval,
    conf.int = cint,
    estimate = c(estimate1,estimate2),
    null.value = null.gam,
    alternative = alternative,
    method = method,
    data.name = dname
  )

  attr(rval, "class") <-"htest"

  return(rval)
} #end of the function
#'
