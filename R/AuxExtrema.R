#AuxExtrema.R;
#Contains the ancillary functions for finding the extrema usually in a restricted fashion
#################################################################

#' @title The furthest points from vertices in each \eqn{CC}-vertex region in a standard basic triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the furthest data points among the data set, \code{Xp}, in each \eqn{CC}-vertex region from the
#' corresponding vertex in the standard basic triangle \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2))}.
#'
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' \code{ch.all.intri} is for checking whether all data points are inside \eqn{T_b} (default is \code{FALSE}).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle.
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}
#' @param ch.all.intri A logical argument for checking whether all data points are inside \eqn{T_b}
#' (default is \code{FALSE}).
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances from furthest points to ..."}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, furthest points from vertices in each vertex region.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \eqn{T_b}.}
#' \item{cent}{The center point used for construction of edge regions.}
#' \item{ncent}{Name of the center, \code{cent}, it is circumcenter \code{"CC"} for this function.}
#' \item{regions}{Vertex regions inside the triangle, \eqn{T_b}, provided as a \code{list}.}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \eqn{T_b}.}
#' \item{dist2ref}{Distances from furthest points in each vertex region to the corresponding vertex.}
#'
#' @seealso \code{\link{fr2vertsCC_VR}}, \code{\link{fr2edgesCM_ERstd.tri}}, \code{\link{kfr2vertsCC_VRbas.tri}} and \code{\link{kfr2vertsCC_VR}}
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
#' n<-20
#'
#' set.seed(1)
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' Ext<-fr2vertsCC_VRbas.tri(Xp,c1,c2)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' f2v<-Ext
#'
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",
#' main="Furthest Points in CC-Vertex Regions \n from the Vertices",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp)
#' points(rbind(f2v$ext),pch=4,col=2)
#'
#' txt<-rbind(Tb,CC,D1,D2,D3)
#' xc<-txt[,1]+c(-.03,.03,0.02,.07,.06,-.05,.01)
#' yc<-txt[,2]+c(.02,.02,.03,.01,.02,.02,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
fr2vertsCC_VRbas.tri <- function(Xp,c1,c2,ch.all.intri=FALSE)
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

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); tri<-rbind(y1,y2,y3)
  CC<-circ.cent.tri(tri)
  D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
  Ds<-rbind(D1,D2,D3)
  L<-matrix(rep(CC,3),ncol=3,byrow=TRUE); R<-Ds

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,tri,boundary=TRUE))
    {stop('not all points in the data set are in the standard basic triangle')}
  }

  n<-nrow(Xp)
  mdt<-rep(0,3); U<-matrix(NA,nrow=3,ncol=2);
  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],tri,boundary = TRUE)$in.tri)
    {
      rv<-rv.bas.triCC(Xp[i,],c1,c2)$rv;
      if (rv==1)
      {d1<-Dist(Xp[i,],y1);
      if (d1>=mdt[1]) {mdt[1]<-d1; U[1,]<-Xp[i,]};
      } else {
        if (rv==2)
        {d2<-Dist(Xp[i,],y2);
        if (d2>=mdt[2]) {mdt[2]<-d2; U[2,]<-Xp[i,]}
        } else {
          d3<-Dist(Xp[i,],y3);
          if (d3>=mdt[3]) {mdt[3]<-d3; U[3,]<-Xp[i,]}
        }}
    }
  }
  row.names(tri)<-c("A","B","C")  #vertex labeling
  Cvec= paste(round(y3,2), collapse=",")
  typ<-paste("Furthest Points in CC-Vertex Regions of the Standard Basic Triangle with Vertices A=(0,0), B=(1,0), and C=(",Cvec,") from its Vertices",sep="")
  description<-"Furthest Points in CC-Vertex Regions of the Standard Basic Triangle from its Vertices \n (Row i corresponds to vertex i for i=1,2,3)"
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the vertices and the furthest points to vertices in each vertex region\n (i-th entry corresponds to vertex i for i=1,2,3)"
  main.txt<-"Furthest Points in CC-Vertex Regions \n from the Vertices"

  Dis<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),ifelse(!is.na(U[3,1]),mdt[3],NA))
  #distances of the furthest points to the vertices
  Regs<-list(vr1=rbind(y1,D3,CC,D2), #regions inside the triangles
             vr2=rbind(y2,D1,CC,D3),
             vr3=rbind(y3,D2,CC,D1)
  )

  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

  supp.type = "Standard Basic Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, #furthest points from vertices in each vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (standard basic triangle here)
    cent=CC, cent.name="CC", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of furthest points to vertices in each vertex region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The furthest points in a data set from vertices in each \eqn{CC}-vertex region in a triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the furthest data points among the data set, \code{Xp}, in each \eqn{CC}-vertex region from the vertex in the
#' triangle, \code{tri}\eqn{=T(A,B,C)}. Vertex region labels/numbers correspond to the row number of the vertex in \code{tri}.
#' \code{ch.all.intri} is for checking whether all data points are inside \code{tri} (default is \code{FALSE}).
#'
#' If some of the data points are not inside \code{tri} and \code{ch.all.intri=TRUE}, then the function yields
#' an error message.
#' If some of the data points are not inside \code{tri} and \code{ch.all.intri=FALSE}, then the function yields
#' the closest points to edges among the data points inside \code{tri} (yields \code{NA} if there are no data points
#' inside \code{tri}).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param ch.all.intri A logical argument (default=\code{FALSE}) to check whether all data points are inside
#' the triangle \code{tri}. So, if it is \code{TRUE},
#' the function checks if all data points are inside the closure of the triangle (i.e., interior and boundary
#' combined) else it does not.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances from furthest points to ..."}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, furthest points from vertices in each \eqn{CC}-vertex region in the
#' triangle \code{tri}.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is the triangle \code{tri} for this function.}
#' \item{cent}{The center point used for construction of edge regions.}
#' \item{ncent}{Name of the center, \code{cent}, it is circumcenter \code{"CC"} for this function}
#' \item{regions}{CC-Vertex regions inside the triangle, \code{tri}, provided as a \code{list}}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{tri}}
#' \item{dist2ref}{Distances from furthest points in each vertex region to the corresponding vertex}
#'
#' @seealso \code{\link{fr2vertsCC_VRbas.tri}}, \code{\link{fr2edgesCM_ERstd.tri}}, \code{\link{kfr2vertsCC_VRbas.tri}} and \code{\link{kfr2vertsCC_VR}}
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
#' Ext<-fr2vertsCC_VR(Xp,Tr)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' f2v<-Ext
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,xlab="",asp=1,ylab="",pch=".",
#' main="Furthest Points in CC-Vertex Regions \n from the Vertices",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp)
#' points(rbind(f2v$ext),pch=4,col=2)
#'
#' txt<-rbind(Tr,CC,Ds)
#' xc<-txt[,1]+c(-.06,.08,.05,.12,-.1,-.1,-.09)
#' yc<-txt[,2]+c(.02,-.02,.05,.0,.02,.06,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' Xp2<-rbind(Xp,c(.2,.4))
#' fr2vertsCC_VR(Xp2,Tr,ch.all.intri = FALSE)
#' #gives an error message if ch.all.intri = TRUE
#' #since not all points in the data set are in the triangle
#' }
#'
#' @export fr2vertsCC_VR
fr2vertsCC_VR <- function(Xp,tri,ch.all.intri=FALSE)
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

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];
  CC<-circ.cent.tri(tri)
  D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
  Ds<-rbind(D1,D2,D3)
  L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
  mdt<-rep(0,3); U<-matrix(NA,nrow=3,ncol=2);

  n<-nrow(Xp)

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,tri,boundary=TRUE))
    {stop('not all points in the data set are in the triangle')}
  }

  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],tri,boundary = TRUE)$in.tri)
    {rv<-rv.triCC(Xp[i,],tri)$rv;
    if (rv==1)
    {d1<-Dist(Xp[i,],y1);
    if (d1>=mdt[1]) {mdt[1]<-d1; U[1,]<-Xp[i,]};
    } else {
      if (rv==2)
      {d2<-Dist(Xp[i,],y2);
      if (d2>=mdt[2]) {mdt[2]<-d2; U[2,]<-Xp[i,]}
      } else {
        d3<-Dist(Xp[i,],y3);
        if (d3>=mdt[3]) {mdt[3]<-d3; U[3,]<-Xp[i,]}
      }}
    }
  }

  row.names(tri)<-c("A","B","C")  #vertex labeling

  Avec= paste(round(y1,2), collapse=","); Bvec= paste(round(y2,2), collapse=","); Cvec= paste(round(y3,2), collapse=",");
  typ<-paste("Furthest Points in CC-Vertex Regions of the Triangle with Vertices A=(",Avec,"), B=(",Bvec,"), and C=(",Cvec,") from its Vertices",sep="")
  description<-"Furthest Points in CC-Vertex Regions of the Triangle from its Vertices \n (Row i corresponds to vertex i for i=1,2,3)"
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the vertices and the furthest points in the vertex regions \n (i-th entry corresponds to vertex i for i=1,2,3)"
  main.txt<-"Furthest Points in CC-Vertex Regions \n from the Vertices"

  Dis<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),ifelse(!is.na(U[3,1]),mdt[3],NA))
  #distances of the furthest points to the vertices
  Regs<-list(vr1=rbind(y1,D3,CC,D2), #regions inside the triangles
             vr2=rbind(y2,D1,CC,D3),
             vr3=rbind(y3,D2,CC,D1)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

  supp.type = "Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, #furthest points from vertices in each vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (triangle here)
    cent=CC, cent.name="CC", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of furthest points to vertices in each vertex region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The \code{k} furthest points from vertices in each \eqn{CC}-vertex region in a standard basic triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the \code{k} furthest data points among the data set, \code{Xp}, in each \eqn{CC}-vertex region from the vertex in the
#' standard basic triangle \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2))}.
#'
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' \code{ch.all.intri} is for checking whether all data points are inside \eqn{T_b} (default is \code{FALSE}).
#' In the extrema, \eqn{ext}, in the output, the first \code{k} entries are the \code{k} furthest points from vertex 1,
#' second \code{k} entries are \code{k} furthest points are from vertex 2, and
#' last \code{k} entries are the \code{k} furthest points from vertex 3
#'  If data size does not allow, \code{NA}'s are inserted for some or all of the \code{k} furthest points for each vertex.
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle.
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}
#' @param k A positive integer. \code{k} furthest data points in each \eqn{CC}-vertex region are to be found if exists, else
#' \code{NA} are provided for (some of) the \code{k} furthest points.
#' @param ch.all.intri A logical argument for checking whether all data points are inside \eqn{T_b}
#' (default is \code{FALSE}).
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A shorter description of the distances as \code{"Distances of k furthest points in the vertex regions to Vertices"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, \code{k} furthest points from vertices in each vertex region.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \eqn{T_b}.}
#' \item{cent}{The center point used for construction of edge regions.}
#' \item{ncent}{Name of the center, \code{cent}, it is circumcenter \code{"CC"} for this function.}
#' \item{regions}{Vertex regions inside the triangle, \eqn{T_b}, provided as a \code{list}.}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \eqn{T_b}.}
#' \item{dist2ref}{Distances from \code{k} furthest points in each vertex region to the corresponding vertex
#' (each row representing a vertex).}
#'
#' @seealso \code{\link{fr2vertsCC_VRbas.tri}}, \code{\link{fr2vertsCC_VR}}, \code{\link{fr2edgesCM_ERstd.tri}}, and \code{\link{kfr2vertsCC_VR}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4; c2<-.6;
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C)
#' n<-20
#' k<-3
#'
#' set.seed(1)
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' Ext<-kfr2vertsCC_VRbas.tri(Xp,c1,c2,k)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' kf2v<-Ext
#'
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",
#' main=paste(k," Furthest Points in CC-Vertex Regions \n from the Vertices",sep=""),
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp)
#' points(kf2v$ext,pch=4,col=2)
#'
#' txt<-rbind(Tb,CC,Ds)
#' xc<-txt[,1]+c(-.03,.03,.02,.07,.06,-.05,.01)
#' yc<-txt[,2]+c(.02,.02,.03,-.02,.02,.03,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
kfr2vertsCC_VRbas.tri <- function(Xp,c1,c2,k,ch.all.intri=FALSE)
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

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)
  CC<-circ.cent.tri(Tb)
  D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
  Ds<-rbind(D1,D2,D3)
  L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,Tb,boundary=TRUE))
    {stop('not all points in the data set are in the standard basic triangle')}
  }

  U1<-U2<-U3<-matrix(NA,nrow=k,ncol=2)
  Dis1<-Dis2<-Dis3<-rep(NA,k)

  n<-nrow(Xp);
  rv<-rep(0,n);
  for (i in 1:n)
  {rv[i]<-rv.bas.triCC(Xp[i,],c1,c2)$rv};

  Xp1<-matrix(Xp[rv==1,],ncol=2);
  Xp2<-matrix(Xp[rv==2,],ncol=2);
  Xp3<-matrix(Xp[rv==3,],ncol=2);

  n1<-nrow(Xp1); n2<-nrow(Xp2); n3<-nrow(Xp3);

  if (n1>0)
  {
    Dis1<-rep(0,n1)
    for (i in 1:n1)
    { Dis1[i]<-Dist(Xp1[i,],y1) }

    ord1<-order(-Dis1)
    K1<-min(k,n1)
    U1[1:K1,]<-Xp1[ord1[1:K1],]
  }

  if (n2>0)
  {
    Dis2<-rep(0,n2)
    for (i in 1:n2)
    { Dis2[i]<-Dist(Xp2[i,],y2) }

    ord2<-order(-Dis2)
    K2<-min(k,n2)
    U2[1:K2,]<-Xp2[ord2[1:K2],]
  }

  if (n3>0)
  {
    Dis3<-rep(0,n3)
    for (i in 1:n3)
    { Dis3[i]<-Dist(Xp3[i,],y3) }
    ord3<-order(-Dis3)
    K3<-min(k,n3)
    U3[1:K3,]<-Xp3[ord3[1:K3],]
  }
  U<-rbind(U1,U2,U3)

  row.names(Tb)<-c("A","B","C")  #vertex labeling
  rn1<-rn2<-rn3<-vector()
  for (i in 1:k) {rn1<-c(rn1,paste(i,". furthest from vertex 1",sep=""));
  rn2<-c(rn2,paste(i,". furthest from vertex 2",sep=""));
  rn3<-c(rn3,paste(i,". furthest from vertex 3",sep=""))}
  row.names(U)<-c(rn1,rn2,rn3)  #extrema labeling

  Cvec= paste(round(y3,2), collapse=",")
  typ<-paste(k, " Furthest Points in CC-Vertex Regions of the Standard Basic Triangle with Vertices A=(0,0), B=(1,0), and C=(",Cvec,") from its Vertices",sep="")
  description<-paste(k, " Furthest Points in CC-Vertex Regions of the Standard Basic Triangle from its Vertices \n (Row i corresponds to vertex i for i=1,2,3)",sep="")
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-paste("Distances between the vertices and the ",k," furthest points in the vertex regions \n (i-th entry corresponds to vertex i for i=1,2,3)",sep="")
  main.txt<-paste(k, " Furthest Points in CC-Vertex Regions \n from the Vertices",sep="")

  Dis<-c(Dis1[1:k],Dis2[1:k],Dis3[1:k])
  #distances of the furthest points to the vertices
  Regs<-list(vr1=rbind(y1,D3,CC,D2), #regions inside the triangles
             vr2=rbind(y2,D1,CC,D3),
             vr3=rbind(y3,D2,CC,D1)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

  supp.type = "Standard Basic Triangle" #name of the support


  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, #k furthest points from vertices in each vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=Tb, supp.type = supp.type, # region of interest for X points, and its type (standard basic triangle here)
    cent=CC, cent.name="CC", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of k furthest points to vertices in each vertex region (earh row corresponds to a vertex)
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The \code{k} furthest points in a data set from vertices in each \eqn{CC}-vertex region in a triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the \code{k} furthest data points among the data set, \code{Xp}, in each \eqn{CC}-vertex region from the vertex in the
#' triangle, \code{tri}\eqn{=T(A,B,C)}, vertices are stacked row-wise. Vertex region labels/numbers correspond to the
#' row number of the vertex in \code{tri}.
#'
#' \code{ch.all.intri} is for checking whether all data points are inside \code{tri} (default is \code{FALSE}).
#' If some of the data points are not inside \code{tri} and \code{ch.all.intri=TRUE}, then the function yields
#' an error message.
#' If some of the data points are not inside \code{tri} and \code{ch.all.intri=FALSE}, then the function yields
#' the closest points to edges among the data points inside \code{tri} (yields \code{NA} if there are no data points
#' inside \code{tri}).
#'
#' In the extrema, \eqn{ext}, in the output, the first \code{k} entries are the \code{k} furthest points from vertex 1,
#' second \code{k} entries are \code{k} furthest points are from vertex 2, and
#' last \code{k} entries are the \code{k} furthest points from vertex 3.
#' If data size does not allow, \code{NA}'s are inserted for some or all of the \code{k} furthest points for each vertex.
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param k A positive integer. \code{k} furthest data points in each \eqn{CC}-vertex region are to be found if exists, else
#' \code{NA} are provided for (some of) the \code{k} furthest points.
#' @param ch.all.intri A logical argument (default=\code{FALSE}) to check whether all data points are inside
#' the triangle \code{tri}. So, if it is \code{TRUE},
#' the function checks if all data points are inside the closure of the triangle (i.e., interior and boundary
#' combined) else it does not.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A shorter description of the distances as \code{"Distances of k furthest points in the vertex regions
#'  to Vertices"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, \code{k} furthest points from vertices in each \eqn{CC}-vertex region in
#' the triangle \code{tri}.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, it is \code{tri} for this function.}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is circumcenter \code{"CC"} for this function.}
#' \item{regions}{Vertex regions inside the triangle, \code{tri}, provided as a \code{list}}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \eqn{T_b}.}
#' \item{dist2ref}{Distances from \code{k} furthest points in each vertex region to the corresponding vertex
#' (each row representing a vertex in \code{tri}). Among the distances the first \code{k} entries are the distances
#' from the \code{k} furthest points from vertex 1 to vertex 1, second \code{k} entries are distances from the \code{k} furthest
#' points from vertex 2 to vertex 2, and the last \code{k} entries are the distances from the \code{k} furthest points
#' from vertex 3 to vertex 3.}
#'
#' @seealso \code{\link{kfr2vertsCC_VRbas.tri}}, \code{\link{fr2vertsCC_VRbas.tri}}, \code{\link{fr2vertsCC_VR}}, and \code{\link{fr2edgesCM_ERstd.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' n<-10  #try also n<-20
#' k<-3
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' Ext<-kfr2vertsCC_VR(Xp,Tr,k)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' Xp2<-rbind(Xp,c(.2,.4))
#' kfr2vertsCC_VR(Xp2,Tr,k)  #try also kfr2vertsCC_VR(Xp2,Tr,k,ch.all.intri = TRUE)
#'
#' kf2v<-Ext
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",
#' main=paste(k," Furthest Points in CC-Vertex Regions \n from the Vertices",sep=""),
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp)
#' points(kf2v$ext,pch=4,col=2)
#'
#' txt<-rbind(Tr,CC,Ds)
#' xc<-txt[,1]+c(-.06,.08,.05,.12,-.1,-.1,-.09)
#' yc<-txt[,2]+c(.02,-.02,.04,.0,.02,.06,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export kfr2vertsCC_VR
kfr2vertsCC_VR <- function(Xp,tri,k,ch.all.intri=FALSE)
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

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];
  CC<-circ.cent.tri(tri)
  D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
  Ds<-rbind(D1,D2,D3)
  L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,tri,boundary=TRUE))
    {stop('not all points in the data set are in the triangle')}
  }

  U1<-U2<-U3<-matrix(NA,nrow=k,ncol=2)
  Dis1<-Dis2<-Dis3<-rep(NA,k)

  n<-nrow(Xp)
  rv<-rep(0,n);

  for (i in 1:n)
  {rv[i]<-rv.triCC(Xp[i,],tri)$rv};

  Xp1<-matrix(Xp[rv==1,],ncol=2);
  Xp2<-matrix(Xp[rv==2,],ncol=2);
  Xp3<-matrix(Xp[rv==3,],ncol=2);

  n1<-nrow(Xp1); n2<-nrow(Xp2); n3<-nrow(Xp3);

  if (n1>0)
  {
    Dis1<-rep(0,n1)
    for (i in 1:n1)
    { Dis1[i]<-Dist(Xp1[i,],y1) }

    ord1<-order(-Dis1)
    K1<-min(k,n1)
    U1[1:K1,]<-Xp1[ord1[1:K1],]
  }

  if (n2>0)
  {
    Dis2<-rep(0,n2)
    for (i in 1:n2)
    { Dis2[i]<-Dist(Xp2[i,],y2) }

    ord2<-order(-Dis2)
    K2<-min(k,n2)
    U2[1:K2,]<-Xp2[ord2[1:K2],]
  }

  if (n3>0)
  {
    Dis3<-rep(0,n3)
    for (i in 1:n3)
    { Dis3[i]<-Dist(Xp3[i,],y3) }
    ord3<-order(-Dis3)
    K3<-min(k,n3)
    U3[1:K3,]<-Xp3[ord3[1:K3],]
  }
  U<-rbind(U1,U2,U3)

  row.names(tri)<-c("A","B","C")  #vertex labeling
  rn1<-rn2<-rn3<-vector()
  for (i in 1:k) {rn1<-c(rn1,paste(i,". furthest from vertex 1",sep=""));
  rn2<-c(rn2,paste(i,". furthest from vertex 2",sep=""));
  rn3<-c(rn3,paste(i,". furthest from vertex 3",sep=""))}
  row.names(U)<-c(rn1,rn2,rn3)  #extrema labeling

  Avec= paste(round(y1,2), collapse=","); Bvec= paste(round(y2,2), collapse=","); Cvec= paste(round(y3,2), collapse=",");
  typ<-paste(k," Furthest Points in CC-Vertex Regions of the Triangle with Vertices A=(",Avec,"), B=(",Bvec,"), and C=(",Cvec,") from its Vertices",sep="")
  description<-paste(k, " Furthest Points in CC-Vertex Regions of the Triangle from its Vertices \n (Row i corresponds to vertex i for i=1,2,3)",sep="")
  txt1<-paste("Vertex labels are A=1, B=2, and C=3 (where vertex i corresponds to row numbers ", k,"(i-1) to ",k,"i in Extremum Points)",sep="")
  txt2<-paste("Distances between the vertices and the ",k," furthest points in the vertex regions \n (i-th entry corresponds to vertex i for i=1,2,3)",sep="")
  main.txt<-paste(k, " Furthest Points in CC-Vertex Regions \n from the Vertices",sep="")

  Dis<-rbind(Dis1[1:k],Dis2[1:k],Dis3[1:k])
  #distances of the furthest points to the vertices
  Regs<-list(vr1=rbind(y1,D3,CC,D2), #regions inside the triangles
             vr2=rbind(y2,D1,CC,D3),
             vr3=rbind(y3,D2,CC,D1)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

  supp.type = "Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, #k furthest points from vertices in each vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (triangle here)
    cent=CC, cent.name="CC", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of k furthest points to vertices in each vertex region (each row corresponds to a vertex)
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points to circumcenter in each \eqn{CC}-vertex region in a standard basic triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to circumcenter, \eqn{CC}, in each \eqn{CC}-vertex region
#' in the standard basic triangle \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2))=}(vertex 1,vertex 2,vertex 3).
#' \code{ch.all.intri} is for checking whether all data points are inside \eqn{T_b} (default is \code{FALSE}).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle.
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}
#' @param ch.all.intri A logical argument for checking whether all data points are inside \eqn{T_b}
#' (default is \code{FALSE}).
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances from closest points to ..."}.}
#' \item{type}{Type of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to \eqn{CC} in each vertex region.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \eqn{T_b}.}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is \code{"CC"} for this function.}
#' \item{regions}{Vertex regions inside the triangle, \eqn{T_b}, provided as a \code{list}.}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \eqn{T_b}.}
#' \item{dist2ref}{Distances from closest points in each vertex region to CC.}
#'
#' @seealso \code{\link{cl2CC.VR}}, \code{\link{cl2edgesVRbas.tri}}, \code{\link{cl2edgesM_VR}}, \code{\link{cl2edgesCM_VR}}
#' and \code{\link{fr2edgesCM_ERstd.tri}}
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
#' n<-15
#'
#' set.seed(1)
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' Ext<-cl2CC.VRbas.tri(Xp,c1,c2)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' c2CC<-Ext
#'
#' CC<-circ.cent.bas.tri(c1,c2)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",
#' main="Closest Points in CC-Vertex Regions \n to the Circumcenter",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp)
#' points(c2CC$ext,pch=4,col=2)
#'
#' txt<-rbind(Tb,CC,Ds)
#' xc<-txt[,1]+c(-.03,.03,.02,.07,.06,-.05,.01)
#' yc<-txt[,2]+c(.02,.02,.03,-.01,.03,.03,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' Xp2<-rbind(Xp,c(.2,.4))
#' cl2CC.VRbas.tri(Xp2,c1,c2,ch.all.intri = FALSE)
#' #gives an error message if ch.all.intri = TRUE
#' #since not all points are in the standard basic triangle
#' }
#'
#' @export
cl2CC.VRbas.tri <- function(Xp,c1,c2,ch.all.intri=FALSE)
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

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)
  CC<-circ.cent.bas.tri(c1,c2)
  D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
  Ds<-rbind(D1,D2,D3)
  L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,Tb,boundary=TRUE))
    {stop('not all points in the data set are in the standard basic triangle')}
  }

  mdt<-c(Dist(y1,CC),Dist(y2,CC),Dist(y3,CC)); #distances from the vertices to CC
  U<-matrix(NA,nrow=3,ncol=2);

  n<-nrow(Xp)
  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],Tb,boundary = TRUE)$in.tri)
    {
      rv<-rv.bas.triCC(Xp[i,],c1,c2)$rv;
      if (rv==1)
      {d1<-Dist(Xp[i,],CC);
      if (d1<=mdt[1]) {mdt[1]<-d1; U[1,]<-Xp[i,]};
      } else {
        if (rv==2)
        {d2<-Dist(Xp[i,],CC);
        if (d2<=mdt[2]) {mdt[2]<-d2; U[2,]<-Xp[i,]}
        } else {
          d3<-Dist(Xp[i,],CC);
          if (d3<=mdt[3]) {mdt[3]<-d3; U[3,]<-Xp[i,]}
        }}
    }
  }

  row.names(Tb)<-c("A","B","C")  #vertex labeling
  Cvec= paste(round(y3,2), collapse=",")
  typ<-paste("Closest Points in CC-Vertex Regions of the Standard Basic Triangle with Vertices A=(0,0), B=(1,0), and C=(",Cvec,") to its Circumcenter",sep="")
  description<-"Closest Points in CC-Vertex Regions of the Standard Basic Triangle to its Circumcenter \n (Row i corresponds to vertex i for i=1,2,3)"
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the Circumcenter and the Closest Points to the Circumcenter in CC-Vertex Regions \n (i-th entry corresponds to vertex i for i=1,2,3)"
  main.txt<-"Closest Points in CC-Vertex Regions \n to the Circumcenter"

  Dis<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),ifelse(!is.na(U[3,1]),mdt[3],NA))
  #distances of the closest points to CC
  Regs<-list(vr1=rbind(y1,D3,CC,D2), #regions inside the triangles
             vr2=rbind(y2,D1,CC,D3),
             vr3=rbind(y3,D2,CC,D1)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

  supp.type = "Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, #closest points to CC in each vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=Tb, supp.type = supp.type, # region of interest for X points, and its type (standard basic triangle here)
    cent=CC, cent.name="CC", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points to CC in each vertex region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points to circumcenter in each \eqn{CC}-vertex region in a triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to circumcenter, \eqn{CC}, in each \eqn{CC}-vertex region
#' in the triangle \code{tri} \eqn{=T(A,B,C)=}(vertex 1,vertex 2,vertex 3).
#'
#' \code{ch.all.intri} is for checking whether all data points are inside \code{tri} (default is \code{FALSE}).
#' If some of the data points are not inside \code{tri} and \code{ch.all.intri=TRUE}, then the function yields
#' an error message.
#' If some of the data points are not inside \code{tri} and \code{ch.all.intri=FALSE}, then the function yields
#' the closest points to \eqn{CC} among the data points in each \eqn{CC}-vertex region of \code{tri} (yields \code{NA} if
#' there are no data points inside \code{tri}).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:mcap2012;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param ch.all.intri A logical argument (default=\code{FALSE}) to check whether all data points are inside
#' the triangle \code{tri}. So, if it is \code{TRUE},
#' the function checks if all data points are inside the closure of the triangle (i.e., interior and boundary
#' combined) else it does not.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances from closest points to CC ..."}}
#' \item{type}{Type of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to \eqn{CC} in each \eqn{CC}-vertex region}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \code{tri}}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is \code{"CC"} for this function}
#' \item{regions}{Vertex regions inside the triangle, \code{tri}, provided as a \code{list}}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{tri}}
#' \item{dist2ref}{Distances from closest points in each \eqn{CC}-vertex region to CC.}
#'
#' @seealso \code{\link{cl2CC.VRbas.tri}}, \code{\link{cl2edgesVRbas.tri}}, \code{\link{cl2edgesM_VR}},
#' \code{\link{cl2edgesCM_VR}} and  \code{\link{fr2edgesCM_ERstd.tri}}
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
#' Ext<-cl2CC.VR(Xp,Tr)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' c2CC<-Ext
#'
#' CC<-circ.cent.tri(Tr)  #the circumcenter
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",asp=1,xlab="",ylab="",
#' main="Closest Points in CC-Vertex Regions \n to the Circumcenter",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(c2CC$ext,pch=4,col=2)
#'
#' txt<-rbind(Tr,CC,Ds)
#' xc<-txt[,1]+c(-.07,.08,.06,.12,-.1,-.1,-.09)
#' yc<-txt[,2]+c(.02,-.02,.03,.0,.02,.06,-.04)
#' txt.str<-c("A","B","C","CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#'
#' Xp2<-rbind(Xp,c(.2,.4))
#' cl2CC.VR(Xp2,Tr,ch.all.intri = FALSE)
#' #gives an error message if ch.all.intri = TRUE since not all points are in the triangle
#' }
#'
#' @export
cl2CC.VR <- function(Xp,tri,ch.all.intri=FALSE)
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

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];
  CC<-circ.cent.tri(tri)
  D1<-(y2+y3)/2; D2<-(y1+y3)/2; D3<-(y1+y2)/2;
  Ds<-rbind(D1,D2,D3)
  L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,tri,boundary=TRUE))
    {stop('not all points in the data set are in the triangle')}
  }

  mdt<-c(Dist(y1,CC),Dist(y2,CC),Dist(y3,CC)); #distances from the vertices to CC
  U<-matrix(NA,nrow=3,ncol=2);

  n<-nrow(Xp)
  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],tri,boundary = TRUE)$in.tri)
    {rv<-rv.triCC(Xp[i,],tri)$rv;
    if (rv==1)
    {d1<-Dist(Xp[i,],CC);
    if (d1<=mdt[1]) {mdt[1]<-d1; U[1,]<-Xp[i,]};
    } else {
      if (rv==2)
      {d2<-Dist(Xp[i,],CC);
      if (d2<=mdt[2]) {mdt[2]<-d2; U[2,]<-Xp[i,]}
      } else {
        d3<-Dist(Xp[i,],CC);
        if (d3<=mdt[3]) {mdt[3]<-d3; U[3,]<-Xp[i,]}
      }}
    }
  }

  Mdt<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),ifelse(!is.na(U[3,1]),mdt[3],NA))

  row.names(tri)<-c("A","B","C")  #vertex labeling

  Avec= paste(round(y1,2), collapse=","); Bvec= paste(round(y2,2), collapse=","); Cvec= paste(round(y3,2), collapse=",");
  typ<-paste("Closest Points in CC-Vertex Regions of the Triangle with Vertices A=(",Avec,"), B=(",Bvec,"), and C=(",Cvec,") to its Circumcenter",sep="")
  description<-"Closest Points in CC-Vertex Regions of the Triangle to its Circumcenter \n (Row i corresponds to vertex i for i=1,2,3)"
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the Circumcenter and the Closest Points to the Circumcenter in CC-Vertex Regions \n (i-th entry corresponds to vertex i for i=1,2,3)"
  main.txt<-"Closest Points in CC-Vertex Regions \n to the Circumcenter"

  Dis<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),ifelse(!is.na(U[3,1]),mdt[3],NA))
  #distances of the closest points to CC
  Regs<-list(vr1=rbind(y1,D3,CC,D2), #regions inside the triangles
             vr2=rbind(y2,D1,CC,D3),
             vr3=rbind(y3,D2,CC,D1)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

  supp.type = "Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, #closest points to CC in each vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (triangle here)
    cent=CC, cent.name="CC", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points to CC in each vertex region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points to center in each vertex region in an interval
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, in each \eqn{M_c}-vertex region
#' i.e., finds the closest points from right and left to \eqn{M_c} among points of the 1D data set \code{Xp} which reside in
#' in the interval \code{int}\eqn{=(a,b)}.
#'
#' \eqn{M_c} is based on the centrality parameter \eqn{c \in (0,1)},
#' so that \eqn{100c} \% of the length of interval is to the left of \eqn{M_c} and \eqn{100(1-c)} \% of the length of the interval
#' is to the right of \eqn{M_c}. That is, for the interval \eqn{(a,b)}, \eqn{M_c=a+c(b-a)}. If there are no points from \code{Xp} to
#' the left of \eqn{M_c} in the interval, then it yields \code{NA}, and likewise for the right of \eqn{M_c} in the interval.
#'
#' See also (\insertCite{ceyhan:metrika-2012;textual}{pcds}).
#'
#' @param Xp A set or \code{vector} of 1D points from which closest points to \eqn{M_c} are found
#' in the interval \code{int}.
#' @param int A \code{vector} of two real numbers representing an interval.
#' @param c A positive real number in \eqn{(0,1)} parameterizing the center inside \code{int}\eqn{=(a,b)}.
#' For the interval, \code{int}\eqn{=(a,b)}, the parameterized center is \eqn{M_c=a+c(b-a)}.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex Labels are \eqn{a=1} and \eqn{b=2} for the interval \eqn{(a,b)}.}
#' \item{txt2}{A short description of the distances as \code{"Distances from ..."}}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to \eqn{M_c} in each vertex region}
#' \item{ind.ext}{The data indices of extrema points, \code{ext}.}
#' \item{X}{The input data vector, \code{Xp}.}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \code{int}.}
#' \item{cent}{The (parameterized) center point used for construction of vertex regions.}
#' \item{ncent}{Name of the (parameterized) center, \code{cent}, it is \code{"Mc"} for this function.}
#' \item{regions}{Vertex regions inside the interval, \code{int}, provided as a list.}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{int}.}
#' \item{dist2ref}{Distances from closest points in each vertex region to \eqn{M_c}.}
#'
#' @seealso \code{\link{cl2CC.VRbas.tri}} and \code{\link{cl2CC.VR}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c<-.4
#' a<-0; b<-10; int<-c(a,b)
#'
#' Mc<-centMc(int,c)
#'
#' nx<-10
#' xr<-range(a,b,Mc)
#' xf<-(xr[2]-xr[1])*.5
#'
#' Xp<-runif(nx,a,b)
#'
#' Ext<-cl2Mc.int(Xp,int,c)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' cMc<-Ext
#'
#' Xlim<-range(a,b,Xp)
#' xd<-Xlim[2]-Xlim[1]
#'
#' plot(cbind(a,0),xlab="",pch=".",
#' main=paste("Closest Points in Mc-Vertex Regions \n to the Center Mc = ",Mc,sep=""),
#'   xlim=Xlim+xd*c(-.05,.05))
#'   abline(h=0)
#' abline(v=c(a,b,Mc),col=c(1,1,2),lty=2)
#' points(cbind(Xp,0))
#' points(cbind(c(cMc$ext),0),pch=4,col=2)
#' text(cbind(c(a,b,Mc)-.02*xd,-0.05),c("a","b",expression(M[c])))
#' }
#'
#' @export cl2Mc.int
cl2Mc.int <- function(Xp,int,c)
{
  if (!is.point(Xp,length(Xp)))
  {stop('Xp must be a 1D vector of numerical entries')}

  if (!is.point(int))
  {stop('int must a numeric vector of length 2')}

  if (!is.point(c,1) || c <= 0 || c >= 1)
  {stop('c must be a scalar in (0,1)')}

  y1<-int[1]; y2<-int[2];
  if (y1>=y2)
  {stop('interval is degenerate or void, left end must be smaller than right end')}

  Mc<-y1+c*(y2-y1)

  ind.int = which((Xp>=y1 & Xp<=y2)) #indices of original data Xp in the interval int
  Xp<-Xp[ind.int] #data in the interval int

  indL<-which(Xp<=Mc) #indices of data in the interval to the left of center
  ind.int.left =ind.int[indL] #indices of original data  in the interval and to the left of center
  XpL =Xp[indL] #data in the interval to the left of center

  U<-rep(NA,2)# closest data points to the center in vertex regions
  ind.ext<-rep(NA,2)  # data indices of U points
  if (length(XpL)>0)
  {ext.indL = which(XpL == max(XpL)) #index of closest point in data left to the center
  U[1]<-XpL[ext.indL] #closest data left to the center
  ind.ext[1]=ind.int.left[ext.indL] #original data index to closest from left to the center
  };

  indR<-which(Xp>Mc) #indices of data in the interval to the right of center
  ind.indR =ind.int[indR] #indices of original data  in the interval and to the right of center
  XpR =Xp[indR] #data in the interval to the right of center

  if (length(XpR)>0)
  {ext.indR = which(XpR == min(XpR)) #index of closest point in data right to the center
  U[2]<-XpR[ext.indR] #closest data right to the center
  ind.ext[2]=ind.indR[ext.indR] #original data index to closest from left to the center
  };

  names(int)<-c("a","b")  #vertex labeling

  int.vec= paste(int, collapse=",");
  typ<-paste("Closest Points in Mc-Vertex Regions of the Interval (a,b) = (",int.vec,") to its Center Mc = ",Mc,sep="")
  description<-"Closest Points in Mc-Vertex Regions of the Interval to its Center \n (i-th entry corresponds to vertex i for i=1,2)"
  txt1<-"Vertex Labels are a=1 and b=2 for the interval (a,b)"
  txt2<-"Distances between the Center Mc and the Closest Points to Mc in Mc-Vertex Regions \n (i-th entry corresponds to vertex i for i=1,2)"
  main.txt<-paste("Closest Points in Mc-Vertex Regions \n to its Center Mc = ",Mc,sep="")

  Dis<-c(Mc-U[1],U[2]-Mc)  #distances of the closest points to Mc
  Regs<-list(vr1=c(int[1],Mc), #regions inside the interval
             vr2=c(Mc,int[2]))

  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-c(Reg.Cent,mean(Regs[[i]]))}

  Reg.names<-c("vr=1","vr=2")  #regions names

  supp.type = "Interval" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, ind.ext=ind.ext, #ext: closest points to Mc in each vertex region, and ind.ext is their data indices
    X=Xp, num.points=length(Xp), #data points and its size
    ROI=int, supp.type = supp.type, # region of interest for X points, and its type (interval here)
    cent=Mc, cent.name="Mc", #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points to Mc in each vertex region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The furthest points in a data set from edges in each \eqn{CM}-edge region in the standard
#' equilateral triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the furthest data points among the data set, \code{Xp}, in each \eqn{CM}-edge region from the edge in the
#' standard equilateral triangle \eqn{T_e=T(A=(0,0),B=(1,0),C=(1/2,\sqrt{3}/2))}.
#'
#' \code{ch.all.intri} is for checking whether all data points are inside \eqn{T_e} (default is \code{FALSE}).
#'
#' See also (\insertCite{ceyhan:Phd-thesis;textual}{pcds}).
#'
#' @param Xp A set of 2D points, some could be inside and some could be outside standard equilateral triangle
#' \eqn{T_e}.
#' @param ch.all.intri A logical argument used for checking whether all data points are inside \eqn{T_e}
#' (default is \code{FALSE}).
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Edge labels as \eqn{AB=3}, \eqn{BC=1}, and \eqn{AC=2} for \eqn{T_e} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Edges"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, furthest points from edges in each edge region.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \eqn{T_e}.}
#' \item{cent}{The center point used for construction of edge regions.}
#' \item{ncent}{Name of the center, \code{cent}, it is center of mass \code{"CM"} for this function.}
#' \item{regions}{Edge regions inside the triangle, \eqn{T_e}, provided as a list.}
#' \item{region.names}{Names of the edge regions as \code{"er=1"}, \code{"er=2"}, and \code{"er=3"}.}
#' \item{region.centers}{Centers of mass of the edge regions inside \eqn{T_e}.}
#' \item{dist2ref}{Distances from furthest points in each edge region to the corresponding edge.}
#'
#' @seealso \code{\link{fr2vertsCC_VRbas.tri}}, \code{\link{fr2vertsCC_VR}}, \code{\link{kfr2vertsCC_VRbas.tri}},
#' \code{\link{kfr2vertsCC_VR}}, and \code{\link{cl2edges.std.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-20
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Ext<-fr2edgesCM_ERstd.tri(Xp)
#' Ext
#' summary(Ext)
#' plot(Ext,asp=1)
#'
#' ed.far<-Ext
#'
#' Xp2<-rbind(Xp,c(.8,.8))
#' fr2edgesCM_ERstd.tri(Xp2)
#' fr2edgesCM_ERstd.tri(Xp2,ch.all.intri = FALSE) #gives error if ch.all.intri = TRUE
#'
#' A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' CM<-(A+B+C)/3
#' p1<-(A+B)/2
#' p2<-(B+C)/2
#' p3<-(A+C)/2
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",
#' main="Furthest Points in CM-Edge Regions \n of Std Equilateral Triangle from its Edges",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' L<-Te; R<-matrix(rep(CM,3),ncol=2,byrow=TRUE)
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(Xp,xlab="",ylab="")
#' points(ed.far$ext,pty=2,pch=4,col="red")
#'
#' txt<-rbind(Te,CM,p1,p2,p3)
#' xc<-txt[,1]+c(-.03,.03,.03,-.06,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,.02,0,0,0)
#' txt.str<-c("A","B","C","CM","re=2","re=3","re=1")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
fr2edgesCM_ERstd.tri <- function(Xp,ch.all.intri=FALSE)
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

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
  Te<-rbind(A,B,C)
  Cent<-(A+B+C)/3; Cname<-"CM"

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,Te,boundary=TRUE))
    {stop('not all points in the data set are in the standard equilateral triangle Te=T((0,0),(1,0),(1/2,sqrt(3)/2))')}
  }

  n<-nrow(Xp)

  D<-rep(0,3)
  xf<-matrix(NA,nrow=3,ncol=2)
  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],Te,boundary=TRUE)$in.tri)
    {
      re<-reTeCM(Xp[i,])$re
      dis<-c((-.5*Xp[i,2]+.8660254040-.8660254040*Xp[i,1]),(-.5*Xp[i,2]+.8660254040*Xp[i,1]),Xp[i,2])
      if ( dis[re] > D[re])
      {
        D[re]<-dis[re]; xf[re,]<-Xp[i,]
      }
    }
  }

  row.names(Te)<-c("A","B","C")  #vertex labeling

  typ<-paste("Furthest Points in the CM-Edge Regions of the Standard Equilateral Triangle T=(A,B,C) with A=(0,0), B=(1,0), and C=(1/2,sqrt(3)/2) from its Edges",sep="")
  description<-"Furthest Points in the CM-Edge Regions of the Standard Equilateral Triangle from its Edges \n (Row i corresponds to edge i for i=1,2,3)"
  txt1<-"Edge Labels are AB=3, BC=1, and AC=2 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the edges and the furthest points in the edge regions \n (i-th entry corresponds to edge i for i=1,2,3)"
  main.txt<-"Furthest Points in the CM-Edge Regions \n of the Standard Equilateral Triangle from its Edges"

  Dis<-c(ifelse(!is.na(xf[1,1]),D[1],NA),ifelse(!is.na(xf[2,1]),D[2],NA),ifelse(!is.na(xf[3,1]),D[3],NA))
  #distances of the furthest points to the edges in corresponding edge regions

  Regs<-list(r1=rbind(A,B,Cent), #regions inside the triangles
             r2=rbind(B,C,Cent),
             r3=rbind(C,A,Cent)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("er=1","er=2","er=3")  #regions names

  supp.type = "Standard Equilateral Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt,
    ext=xf, #furthest points from edges in each edge region
    X=Xp, num.points=n, #data points and its size
    ROI=Te, supp.type = supp.type, # region of interest for X points, and its type (std eq triangle here)
    cent=Cent, cent.name=Cname, #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of furthest points to edges in each edge region
  )
  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points in a data set to edges in the standard equilateral triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest points from the 2D data set, \code{Xp}, to the edges in the
#' standard equilateral triangle \eqn{T_e=T(A=(0,0),B=(1,0),C=(1/2,\sqrt{3}/2))}.
#'
#' \code{ch.all.intri} is for checking whether all data points are inside \eqn{T_e} (default is \code{FALSE}).
#'
#' If some of the data points are not inside \eqn{T_e} and \code{ch.all.intri=TRUE}, then the function yields
#' an error message.
#' If some of the data points are not inside \eqn{T_e} and \code{ch.all.intri=FALSE}, then the function yields
#' the closest points to edges among the data points inside \eqn{T_e} (yields \code{NA} if there are no data points
#' inside \eqn{T_e}).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:masa-2007;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param ch.all.intri A logical argument (default=\code{FALSE}) to check whether all data points are inside
#' the standard equilateral triangle \eqn{T_e}. So, if it is \code{TRUE},
#' the function checks if all data points are inside the closure of the triangle (i.e., interior and boundary
#' combined) else it does not.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Edge labels as \eqn{AB=3}, \eqn{BC=1}, and \eqn{AC=2} for \eqn{T_e} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Edges ..."}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, i.e., closest points to edges}
#' \item{X}{The input data, \code{Xp}, which can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, i.e., the standard equilateral triangle \eqn{T_e}}
#' \item{cent}{The center point used for construction of edge regions, not required for this extrema,
#' hence it is \code{NULL} for this function}
#' \item{ncent}{Name of the center, \code{cent}, not required for this extrema, hence it is \code{NULL} for this function}
#' \item{regions}{Edge regions inside the triangle, \eqn{T_e}, not required for this extrema, hence it is \code{NULL}
#' for this function}
#' \item{region.names}{Names of the edge regions, not required for this extrema, hence it is \code{NULL}
#' for this function}
#' \item{region.centers}{Centers of mass of the edge regions inside \eqn{T_e}, not required for this extrema,
#' hence it is \code{NULL} for this function}
#' \item{dist2ref}{Distances from closest points in each edge region to the corresponding edge}
#'
#' @seealso \code{\link{cl2edgesVRbas.tri}}, \code{\link{cl2edgesM_VR}}, \code{\link{cl2edgesCM_VR}} and \code{\link{fr2edgesCM_ERstd.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-20  #try also n<-100
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Ext<-cl2edges.std.tri(Xp)
#' Ext
#' summary(Ext)
#' plot(Ext,asp=1)
#'
#' ed.clo<-Ext
#'
#' A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' CM<-(A+B+C)/3
#' p1<-(A+B)/2
#' p2<-(B+C)/2
#' p3<-(A+C)/2
#'
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp,xlab="",ylab="")
#' points(ed.clo$ext,pty=2,pch=4,col="red")
#'
#' txt<-rbind(Te,p1,p2,p3)
#' xc<-txt[,1]+c(-.03,.03,.03,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.02,0,0,0)
#' txt.str<-c("A","B","C","re=1","re=2","re=3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cl2edges.std.tri
cl2edges.std.tri <- function(Xp,ch.all.intri=FALSE)
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

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2);
  Te<-rbind(A,B,C)
  Cent<-c()
  Cname<-NULL

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,Te,boundary=TRUE))
    {stop('not all points in the data set are in the standard equilateral triangle Te=T((0,0),(1,0),(1/2,sqrt(3)/2))')}
  }

  n<-nrow(Xp)

  D<-rep(0.8660254,3); #distance from a vertex to the opposite edge in Te
  xc<-matrix(NA,nrow=3,ncol=2)
  for (i in 1:n)
  {
    if (in.triangle(Xp[i,],Te,boundary=TRUE)$in.tri)
    {
      dis<-c((-.5*Xp[i,2]+.8660254040-.8660254040*Xp[i,1]),(-.5*Xp[i,2]+.8660254040*Xp[i,1]),Xp[i,2])
      for (j in 1:3)
      {
        if (dis[j]<D[j])
        {D[j]<-dis[j]; xc[j,]<-Xp[i,]}
      }
    }
  }

  row.names(Te)<-c("A","B","C")  #vertex labeling
  typ<-"Closest Points in the Standard Equilateral Triangle Te=T(A,B,C) with Vertices A=(0,0), B=(1,0), and C=(1/2,sqrt(3)/2) to its Edges"
  txt1<-"Edge labels are AB=3, BC=1, and AC=2 (correspond to row number in Extremum Points)"
  txt2<-"Distances between Edges and the Closest Points in the Standard Equilateral Triangle \n (Row i corresponds to edge i for i=1,2,3)"
  description<-"Closest Points in the Standard Equilateral Triangle to its Edges  \n (Row i corresponds to edge i for i=1,2,3)  "
  main.txt<-"Closest Points in Standard Equilateral Triangle \n to its Edges"

  Dis<-c(ifelse(!is.na(xc[1,1]),D[1],NA),ifelse(!is.na(xc[2,1]),D[2],NA),ifelse(!is.na(xc[3,1]),D[3],NA))
  #distances of the closest points to the edges in \eqn{T_e}

  Regs<-Reg.Cent<-Reg.names<-c()
  supp.type = "Standard Equilateral Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt,
    ext=xc, #closest points to edges in the std eq triangle
    X=Xp, num.points=n, #data points
    ROI=Te, supp.type = supp.type, # region of interest for X points, and its type (std eq triangle here)
    cent=Cent, cent.name=Cname, #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points to edges in each edge region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points among a data set in the standard equilateral triangle
#' to the median lines in the six  half edge regions
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the six closest points among the data set, \code{Xp}, in the standard equilateral triangle
#' \eqn{T_e=T(A=(0,0),B=(1,0),C=(1/2,\sqrt{3}/2))} in half edge regions. In particular,
#' in regions \eqn{r_1} and \eqn{r_6}, it finds the closest point in each region to the line segment \eqn{[A,CM]}
#' in regions \eqn{r_2} and \eqn{r_3}, it finds the closest point in each region to the line segment \eqn{[B,CM]} and
#' in regions \eqn{r_4} and \eqn{r_5}, it finds the closest point in each region to the line segment \eqn{[C,CM]}
#' where \eqn{CM=(A+B+C)/3} is the center of mass.
#'
#' See the example for this function or example for
#' \code{rel.six.Te} function. If there is no data point in region \eqn{r_i}, then it returns "\code{NA} \code{NA}" for \eqn{i}-th row in the extrema.
#' \code{ch.all.intri} is for checking whether all data points are in \eqn{T_e} (default is \code{FALSE}).
#'
#' @param Xp A set of 2D points among which the closest points in the standard equilateral triangle
#' to the median lines in 6 half edge regions.
#' @param ch.all.intri A logical argument for checking whether all data points are in \eqn{T_e}
#' (default is \code{FALSE}).
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Region labels as r1-r6 (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Line Segments (A,CM), (B,CM), and (C,CM) in the six regions r1-r6"}.}
#' \item{type}{Type of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points in each of regions \code{r1-r6} to the line segments
#' joining vertices to the center of mass, \eqn{CM}.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \eqn{T_e}.}
#' \item{cent}{The center point used for construction of edge regions.}
#' \item{ncent}{Name of the center, \code{cent}, it is center of mass \code{"CM"} for this function.}
#' \item{regions}{The six regions, \code{r1-r6} and edge regions inside the triangle, \eqn{T_e}, provided as a list.}
#' \item{region.names}{Names of the regions as \code{"r1"}-\code{"r6"} and names of the edge regions as \code{"er=1"}, \code{"er=2"}, and \code{"er=3"}.}
#' \item{region.centers}{Centers of mass of the regions \code{r1-r6} and of edge regions inside \eqn{T_e}.}
#' \item{dist2ref}{Distances from closest points in each of regions \code{r1-r6} to the line segments
#' joining vertices to the center of mass, \eqn{CM}.}
#'
#' @seealso \code{\link{rel.six.Te}} and \code{\link{cl2edges.std.tri}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' n<-20  #try also n<-100
#' Xp<-runif.std.tri(n)$gen.points
#'
#' Ext<-six.ext(Xp)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' sixt<-Ext
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
#' Xlim<-range(Te[,1],Xp[,1])
#' Ylim<-range(Te[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' L<-Te; R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' polygon(rbind(h1,h2,h3,h4,h5,h6))
#' points(Xp)
#' points(sixt$ext,pty=2,pch=4,col="red")
#'
#' txt<-rbind(Te,r1,r2,r3,r4,r5,r6)
#' xc<-txt[,1]+c(-.02,.02,.02,0,0,0,0,0,0)
#' yc<-txt[,2]+c(.02,.02,.03,0,0,0,0,0,0)
#' txt.str<-c("A","B","C","1","2","3","4","5","6")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
six.ext <- function(Xp,ch.all.intri=FALSE)
{
  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)
  Cent<-(A+B+C)/3; Cname<-"CM"
  D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
  Ds<-rbind(D1,D2,D3)

  if (ch.all.intri==TRUE)
  {
    if (!in.tri.all(Xp,Te,boundary=TRUE))
    {stop('not all points in the data set are in the standard equilateral triangle Te=T((0,0),(1,0),(1/2,sqrt(3)/2))')}
  }

  if (is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  n<-nrow(Xp)

  D<-rep(0.5773503,6); #distance from CM to each of the vertices in CM
  xc<-matrix(NA,nrow=6,ncol=2)
  for (i in 1:n)
  {rel<-rel.six.Te(Xp[i,])
  if (!is.na(rel))
  {x<-Xp[i,1]; y<-Xp[i,2];
  dis<-c((-0.8660254042*y + 0.5*x),(-.8660254042*y+.5-.5*x),(0.8660254042*y - 0.5+ 0.5*x),(x-.5),(.5-x),
         (.8660254042*y-.5000000003*x))
  if ( dis[rel] < D[rel])
  {
    D[rel]<-dis[rel]; xc[rel,]<-Xp[i,]
  }
  }
  }

  row.names(Te)<-c("A","B","C")  #vertex labeling
  row.names(xc)<-c("closest to line segment (A,CM) in region r1:",
                   "closest to line segment (B,CM) in region r2:",
                   "closest to line segment (B,CM) in region r3:",
                   "closest to line segment (C,CM) in region r4:",
                   "closest to line segment (C,CM) in region r5:",
                   "closest to line segment (A,CM) in region r6:")  #extrema labeling
  typ<-"Closest Points to Line Segments (A,CM), (B,CM), and (C,CM), in the six regions r1-r6 in the Standard Equilateral Triangle with Vertices A=(0,0), B=(1,0), and C=(1/2,sqrt(3)/2) and Center of Mass CM"
  txt1<-"Region labels are r1-r6 (corresponding to row number in Extremum Points)"
  txt2<-"Distances to Line Segments (A,CM), (B,CM), and (C,CM) in the six regions r1-r6"
  description<-"Closest Points to Line Segments (A,CM), (B,CM), and (C,CM) in the Six Regions r1-r6 in the Standard Equilateral Triangle"
  main.txt<-paste("Closest Points to Line Segments \n (A,CM), (B,CM), and (C,CM) in the Regions r1-r6 \n in the Standard Equilateral Triangle")

  h1<-c(1/2,sqrt(3)/18); h2<-c(2/3, sqrt(3)/9); h3<-c(2/3, 2*sqrt(3)/9);
  h4<-c(1/2, 5*sqrt(3)/18); h5<-c(1/3, 2*sqrt(3)/9); h6<-c(1/3, sqrt(3)/9);
  Regs<-list(r1=rbind(h6,h1,Cent), #regions inside the triangles
             r2=rbind(h1,h2,Cent),
             r3=rbind(h2,h3,Cent),
             r4=rbind(h3,h4,Cent),
             r5=rbind(h4,h5,Cent),
             r6=rbind(h5,h6,Cent),
             reg1=rbind(A,B,Cent),
             reg2=rbind(A,C,Cent),
             reg3=rbind(B,C,Cent)
  )

  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("r1","r2","r3","r4","r5","r6"," "," ")  #regions names

  supp.type = "Standard Equilateral Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=xc, #closest points to line segments joining vertices to CM in each of regions r1-r6.
    X=Xp, num.points=n, #data points and its size
    ROI=Te, supp.type = supp.type, # region of interest for X points, and its type (std eq triangle here)
    cent=Cent, cent.name=Cname, #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=D #distances of closest points to line segments joining vertices to CM in each of regions r1-r6.
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points among a data set in the vertex regions to the corresponding edges in a standard basic triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to edge \eqn{i} in \code{M}-vertex region \eqn{i} for \eqn{i=1,2,3}
#' in the standard basic triangle \eqn{T_b=T(A=(0,0),B=(1,0),C=(c_1,c_2))} where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3}, and corresponding edge labels are \eqn{BC=1}, \eqn{AC=2}, and \eqn{AB=3}.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the standard basic triangle \eqn{T_b}
#' or based on the circumcenter of \eqn{T_b}.
#'
#' Any given triangle can be mapped to the standard basic triangle by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling, preserving uniformity of the points in the
#' original triangle. Hence standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param c1,c2 Positive real numbers which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges; \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle \eqn{T_b} or the circumcenter of \eqn{T_b}.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Edges in the Respective \eqn{M}-Vertex Regions"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to edges in the corresponding vertex region.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \eqn{T_b}.}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is \code{"M"} or \code{"CC"} for this function}
#' \item{regions}{Vertex regions inside the triangle, \eqn{T_b}.}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \eqn{T_b}.}
#' \item{dist2ref}{Distances of closest points in the vertex regions to corresponding edges.}
#'
#' @seealso \code{\link{cl2edgesCM_VR}}, \code{\link{cl2edgesM_VR}}, and \code{\link{cl2edges.std.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' c1<-.4;  c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C);
#'
#' set.seed(1)
#' n<-20
#' Xp<-runif.bas.tri(n,c1,c2)$g
#'
#' M<-as.numeric(runif.bas.tri(1,c1,c2)$g)  #try also M<-c(.6,.3)
#'
#' Ext<-cl2edgesVRbas.tri(Xp,c1,c2,M)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' cl2e<-Ext
#'
#' Ds<-cent2edges.bas.tri(c1,c2,M)
#'
#' Xlim<-range(Tb[,1],Xp[,1])
#' Ylim<-range(Tb[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tb,pch=".",xlab="",ylab="",main="Closest Points in M-Vertex Regions \n to the Opposite Edges",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tb)
#' points(Xp,pch=1,col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(cl2e$ext,pch=3,col=2)
#'
#' xc<-Tb[,1]+c(-.02,.02,0.02)
#' yc<-Tb[,2]+c(.02,.02,.02)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(-.02,.04,-.03,0)
#' yc<-txt[,2]+c(-.02,.02,.02,-.03)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export
cl2edgesVRbas.tri <- function(Xp,c1,c2,M)
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

  if (!is.point(c1,1) || !is.point(c2,1))
  {stop('c1 and c2 must be scalars')}

  if (c1<0 || c1>1/2 || c2<=0 || (1-c1)^2+c2^2 >1)
  {stop('c1 must be in [0,1/2], c2 > 0 and (1-c1)^2+c2^2 <= 1')}

  if (!is.point(M) && !is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates')}

  y1<-c(0,0); y2<-c(1,0); y3<-c(c1,c2); Tb<-rbind(y1,y2,y3)

  if (dimension(M)==3)
  {M<-bary2cart(M,Tb)}

  if (isTRUE(all.equal(M,circ.cent.tri(Tb)))==FALSE & in.triangle(M,Tb,boundary=FALSE)$in.tri==FALSE)
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  if (isTRUE(all.equal(M,circ.cent.tri(Tb)))==TRUE)
  {
    res<-cl2edgesCC_VR(Xp,Tb)
    cent.name<-"CC"
  } else
  {
    cent.name<-"M"
    Ds<-cent2edges.bas.tri(c1,c2,M)
    D1<-Ds[1,]; D2<-Ds[2,]; D3<-Ds[3,];
    L<-rbind(M,M,M); R<-Ds

    if (in.triangle(M,Tb,boundary=FALSE)$in.tri==FALSE)
    {stop('center is not in the interior of the standard basic triangle')}

    mdt<-rep(1,3); #maximum distance from a point in the basic tri to its vertices (which is larger than distances to its edges)
    U<-matrix(NA,nrow=3,ncol=2);

    Xp<-matrix(Xp,ncol=2)
    n<-nrow(Xp)
    for (i in 1:n)
    {if (in.triangle(Xp[i,],Tb,boundary = TRUE)$in.tri)
    {rv<-rv.bas.tri.cent(Xp[i,],c1,c2,M)$rv;
    if (rv==1)
    {d1<-dist.point2line(Xp[i,],y2,y3)$dis;
    if (d1<=mdt[1]) {mdt[1]<-d1; U[1,]<-Xp[i,]};
    } else {
      if (rv==2)
      {d2<-dist.point2line(Xp[i,],y1,y3)$dis;
      if (d2<=mdt[2]) {mdt[2]<-d2; U[2,]<-Xp[i,]}
      } else {
        d3<-Xp[i,2];
        if (d3<=mdt[3]) {mdt[3]<-d3; U[3,]<-Xp[i,]}
      }}
    }
    }

    row.names(Tb)<-c("A","B","C")  #vertex labeling
    Cvec= paste(round(y3,2), collapse=",")
    typ<-paste("Closest Points in M-Vertex Regions of the Standard Basic Triangle with Vertices A=(0,0), B=(1,0), and C=(",Cvec,") \n to the Opposite Edges",sep="")
    description<-"Closest Points in M-Vertex Regions of the Standard Basic Triangle to its Edges \n (Row i corresponds to vertex region i for i=1,2,3)"
    txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
    txt2<-"Distances between the Edges and the Closest Points to the Edges in M-Vertex Regions \n (i-th entry corresponds to vertex region i for i=1,2,3)"
    main.txt<-"Closest Points in M-Vertex Regions \n to the Opposite Edges"

    Dis<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),ifelse(!is.na(U[3,1]),mdt[3],NA))
    #distances of the closest points to the edges in the respective vertex regions
    Regs<-list(vr1=rbind(y1,D3,M,D2), #regions inside the triangles
               vr2=rbind(y2,D1,M,D3),
               vr3=rbind(y3,D2,M,D1)
    )
    Reg.Cent<-vector()
    for (i in 1:length(Regs))
    { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

    Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

    supp.type = "Standard Basic Triangle" #name of the support

    res<-list(
      txt1=txt1, txt2=txt2,
      type=typ, desc=description,
      mtitle=main.txt, #main label in the plot
      ext=U, #closest points to edges in each associated vertex region
      X=Xp, num.points=n, #data points and its size
      ROI=Tb, supp.type = supp.type, # region of interest for X points, and its type (standard basic triangle here)
      cent=M, cent.name=cent.name, #center and center name
      regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
      dist2ref=Dis #distances of closest points in vertex regions to the corresponding edges
    )

    class(res)<-"Extrema"
    res$call <-match.call()
  }
  res
} #end of the function
#'

#################################################################

#' @title The closest points among a data set in the vertex regions to the respective edges in a triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to edge \eqn{i} in \code{M}-vertex region \eqn{i} for \eqn{i=1,2,3}
#' in the triangle \code{tri}\eqn{=T(A,B,C)}. Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3}, and corresponding edge labels are
#' \eqn{BC=1}, \eqn{AC=2}, and \eqn{AB=3}.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri}.
#'
#' Two methods of finding these extrema are provided in the function, which can be chosen in the logical argument \code{alt},
#' whose default is \code{alt=FALSE}.
#' When \code{alt=FALSE}, the function sequentially finds the vertex region of the data point and then updates the minimum distance
#' to the opposite edge and the relevant extrema objects, and when \code{alt=TRUE}, it first partitions the data set according which vertex regions they reside, and
#' then finds the minimum distance to the opposite edge and the relevant extrema on each partition.
#' Both options yield equivalent results for the extrema points and indices, with the default being slightly ~ 20% faster.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri} or the circumcenter of \code{tri};
#' which may be entered as "CC" as well;
#' @param alt A logical argument for alternative method of finding the closest points to the edges, default \code{alt=FALSE}.
#' When \code{alt=FALSE}, the function sequentially finds the vertex region of the data point and then the minimum distance
#' to the opposite edge and the relevant extrema objects, and when \code{alt=TRUE}, it first partitions the data set according which vertex regions they reside, and
#' then finds the minimum distance to the opposite edge and the relevant extrema on each partition.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Edges in the Respective \eqn{M}-Vertex Regions"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to edges in the respective vertex region.}
#' \item{ind.ext}{The data indices of extrema points, \code{ext}.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \code{tri}}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is \code{"M"} or \code{"CC"} for this function}
#' \item{regions}{Vertex regions inside the triangle, \code{tri}, provided as a \code{list}}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{tri}}
#' \item{dist2ref}{Distances of closest points in the \code{M}-vertex regions to corresponding edges.}
#'
#' @seealso \code{\link{cl2edgesVRbas.tri}}, \code{\link{cl2edgesCM_VR}}, and \code{\link{cl2edges.std.tri}}
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
#' n<-20  #try also n<-100
#'
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' M<-as.numeric(runif.tri(1,Tr)$g)  #try also M<-c(1.6,1.0)
#'
#' Ext<-cl2edgesM_VR(Xp,Tr,M)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' cl2e<-Ext
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
#' plot(Tr,pch=".",xlab="",ylab="",main="Closest Points in M-Vertex Regions \n to the Opposite Edges",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#' points(Xp,pch=1,col=1)
#' L<-rbind(M,M,M); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(cl2e$ext,pch=3,col=2)
#'
#' xc<-Tr[,1]+c(-.02,.03,.02)
#' yc<-Tr[,2]+c(.02,.02,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' txt<-rbind(M,Ds)
#' xc<-txt[,1]+c(-.02,.05,-.02,-.01)
#' yc<-txt[,2]+c(-.03,.02,.08,-.07)
#' txt.str<-c("M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cl2edgesM_VR
cl2edgesM_VR <- function(Xp,tri,M,alt=FALSE)
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
  {stop('M must be a numeric 2D point for Cartesian coordinates or 3D point for barycentric coordinates
          or the circumcenter "CC" ')}

  CC = circ.cent.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (dimension(M)==3)
  {M<-bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) || in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('center is not the circumcenter or not in the interior of the triangle')}

  y1<-tri[1,]; y2<-tri[2,]; y3<-tri[3,];

  if (isTRUE(all.equal(M,CC))==TRUE)
  {
    res<-cl2edgesCC_VR(Xp,tri)
    cent.name<-"CC"
  } else
  {
    cent.name<-"M"
    Ds<-cent2edges.tri(tri,M)
    D1<-Ds[1,]; D2<-Ds[2,]; D3<-Ds[3,];
    L<-rbind(M,M,M); R<-Ds

    Xp<-matrix(Xp,ncol=2)
    n<-nrow(Xp);

    if (alt)
    {  #A<-tri[1,]; B<-tri[2,]; C<-tri[3,], so A=y1, B=y2, and C=y3

      clBC<-NA; clAC<-NA; clAB<-NA
      VRdt<-rel.verts.tri.cent(Xp,tri,M)$rv
      indA = which(VRdt==1); indB = which(VRdt==2); indC = which(VRdt==3);
      dtA<-matrix(Xp[indA,],ncol=2); dtB<-matrix(Xp[indB,],ncol=2); dtC<-matrix(Xp[indC,],ncol=2)
      distA<-distB<-distC<-vector()
      nA<-nrow(dtA); nB<-nrow(dtB); nC<-nrow(dtC);

      if (nA>0)
      {
        for (i in 1:nA)
        {distA<-c(distA,dist.point2line(dtA[i,],y2,y3)$dis)};
        minA = min(distA)
        cl.ind.dtA = which(distA==minA)
        clBC<-dtA[cl.ind.dtA,]; cl.indA = indA[cl.ind.dtA]
      }

      if (nB>0)
      {
        for (i in 1:nB)
        {distB<-c(distB,dist.point2line(dtB[i,],y1,y3)$dis)};
        minB = min(distB)
        cl.ind.dtB = which(distB==minB)
        clAC<-dtB[cl.ind.dtB,]; cl.indB = indA[cl.ind.dtB]
      }

      if (nC>0)
      {
        for (i in 1:nC)
        {distC<-c(distC,dist.point2line(dtC[i,],y1,y2)$dis)};
        minC = min(distC)
        cl.ind.dtC = which(distC==min(distC))
        clAB<-dtC[cl.ind.dtC,]; cl.indC = indC[cl.ind.dtC]
      }
      U <-rbind(clBC,clAC,clAB); row.names(U)=NULL
      ind.ext = c(cl.indA,cl.indB,cl.indC)
      min.d2e = c(minA,minB,minC) #min distances to edges
    } else
    {
      min.d2e<-c(dist.point2line(y1,y2,y3)$dis,dist.point2line(y2,y1,y3)$dis,dist.point2line(y3,y1,y2)$dis); #distances from each vertex to the opposite edge
      U<-matrix(NA,nrow=3,ncol=2); # closest data points in vertex regions to their respective (opposite to the vertex) edges
      ind.ext<-rep(NA,3)  # data indices of the closest points in vertex regions to their respective (opposite to the vertex) edges

      for (i in 1:n)
      {if (in.triangle(Xp[i,],tri,boundary = TRUE)$in.tri)
      {
        rv<-rv.tri.cent(Xp[i,],tri,M)$rv
        if (rv==1)
        {d1<-dist.point2line(Xp[i,],y2,y3)$dis;
        if (d1<=min.d2e[1]) {min.d2e[1]<-d1; U[1,]<-Xp[i,]; ind.ext[1]=i};
        } else {
          if (rv==2)
          {d2<-dist.point2line(Xp[i,],y1,y3)$dis;
          if (d2<=min.d2e[2]) {min.d2e[2]<-d2; U[2,]<-Xp[i,]; ind.ext[2]=i}
          } else {
            d3<-dist.point2line(Xp[i,],y1,y2)$dis;
            if (d3<=min.d2e[3]) {min.d2e[3]<-d3; U[3,]<-Xp[i,]; ind.ext[3]=i}
          }}
      }
      }
    }

    row.names(tri)<-c("A","B","C")  #vertex labeling

    Avec= paste(round(y1,2), collapse=","); Bvec= paste(round(y2,2), collapse=","); Cvec= paste(round(y3,2), collapse=",");
    typ<-paste("Closest Points in M-Vertex Regions of the Triangle with Vertices A=(",Avec,"), B=(",Bvec,"), and C=(",Cvec,")\n to Opposite Edges",sep="")
    description<-"Closest Points in M-Vertex Regions of the Triangle to its Edges \n (Row i corresponds to vertex region i for i=1,2,3)"
    txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
    txt2<-"Distances between the Edges and the Closest Points to the Edges in M-Vertex Regions \n (i-th entry corresponds to vertex region i for i=1,2,3)"
    main.txt<-"Closest Points in M-Vertex Regions \n to the Opposite Edges"

    Dis<-c(ifelse(!is.na(U[1,1]),min.d2e[1],NA),ifelse(!is.na(U[2,1]),min.d2e[2],NA),ifelse(!is.na(U[3,1]),min.d2e[3],NA))
    #distances of the closest points to the edges in corresponding vertex regions
    Regs<-list(vr1=rbind(y1,D3,M,D2), #regions inside the triangles
               vr2=rbind(y2,D1,M,D3),
               vr3=rbind(y3,D2,M,D1)
    )
    Reg.Cent<-vector()
    for (i in 1:length(Regs))
    { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

    Reg.names<-c("vr=1","vr=2","vr=3")  #regions names

    supp.type = "Triangle" #name of the support

    res<-list(
      txt1=txt1, txt2=txt2,
      type=typ, desc=description,
      mtitle=main.txt, #main label in the plot
      ext=U, ind.ext=ind.ext, #ext: closest points to edges in each associated vertex region, and ind.ext is their data indices
      X=Xp, num.points=n, #data points and its size
      ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (triangle here)
      cent=M, cent.name=cent.name, #center and center name
      regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
      dist2ref=Dis #distances of closest points in vertex regions to the corresponding edges
    )

    class(res)<-"Extrema"
    res$call <-match.call()
  }
  res
} #end of the function
#'

#################################################################

#' @title The closest points in a data set to edges in each \eqn{CM}-vertex region in a triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to edge \eqn{j} in \eqn{CM}-vertex region \eqn{j} for \eqn{j=1,2,3}
#' in the triangle, \code{tri}\eqn{=T(A,B,C)}, where \eqn{CM} stands for center of mass. Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3},
#' and corresponding edge labels are \eqn{BC=1}, \eqn{AC=2}, and \eqn{AB=3}.
#' Function yields \code{NA} if there are no data points in a \eqn{CM}-vertex region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:masa-2007,ceyhan:comp-geo-2010,ceyhan:dom-num-NPE-Spat2011;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Edges in the Respective
#' CM-Vertex Regions"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to edges in the respective vertex region.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \code{tri}}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is \code{"CM"} for this function}
#' \item{regions}{Vertex regions inside the triangle, \code{tri}, provided as a \code{list}}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{tri}}
#' \item{dist2ref}{Distances of closest points in the vertex regions to corresponding edges}
#'
#' @seealso \code{\link{cl2edgesVRbas.tri}}, \code{\link{cl2edgesCC_VR}}, \code{\link{cl2edgesM_VR}},
#' and \code{\link{cl2edges.std.tri}}
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
#' n<-20  #try also n<-100
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' Ext<-cl2edgesCM_VR(Xp,Tr)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' cl2e<-Ext
#'
#' CM<-(A+B+C)/3;
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1])
#' Ylim<-range(Tr[,2],Xp[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,pch=".",xlab="",ylab="",main="Closest Points in CM-Vertex Regions \n to the Opposite Edges",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#'
#' xc<-Tr[,1]+c(-.02,.02,.02)
#' yc<-Tr[,2]+c(.02,.02,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' points(Xp,pch=1,col=1)
#' L<-matrix(rep(CM,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(cl2e$ext,pch=3,col=2)
#'
#' txt<-rbind(CM,Ds)
#' xc<-txt[,1]+c(-.04,.04,-.03,0)
#' yc<-txt[,2]+c(-.05,.04,.06,-.08)
#' txt.str<-c("CM","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cl2edgesCM_VR
cl2edgesCM_VR <- function(Xp,tri)
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

  n<-nrow(Xp)
  A<-tri[1,]; B<-tri[2,]; C<-tri[3,]
  Cent<-(A+B+C)/3; Cname<-"CM"
  D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;

  clBC<-NA; clAC<-NA; clAB<-NA
  VRdt<-rel.verts.triCM(Xp,tri)$rv
  dtA<-matrix(Xp[VRdt==1 & !is.na(VRdt),],ncol=2);
  dtB<-matrix(Xp[VRdt==2 & !is.na(VRdt),],ncol=2);
  dtC<-matrix(Xp[VRdt==3 & !is.na(VRdt),],ncol=2)
  distA<-distB<-distC<-vector()
  nA<-nrow(dtA); nB<-nrow(dtB); nC<-nrow(dtC);

  if (nA>0)
  {
    for (i in 1:nA)
    {distA<-c(distA,dist.point2line(dtA[i,],B,C)$dis)};
    clBC<-dtA[distA==min(distA),]
  }

  if (nB>0)
  {
    for (i in 1:nB)
    {distB<-c(distB,dist.point2line(dtB[i,],A,C)$dis)};
    clAC<-dtB[distB==min(distB),]
  }

  if (nC>0)
  {
    for (i in 1:nC)
    {distC<-c(distC,dist.point2line(dtC[i,],A,B)$dis)};
    clAB<-dtC[distC==min(distC),]
  }
  ce<-rbind(clBC,clAC,clAB)
  row.names(ce)<-c()

  Dis<-c(ifelse(is.numeric(distA),min(distA),NA),ifelse(is.numeric(distB),min(distB),NA),
             ifelse(is.numeric(distC),min(distC),NA))

  row.names(tri)<-c("A","B","C")  #vertex labeling
  txt<-"Edge labels are AB=3, BC=1, and AC=2 (correspond to row number in ce)"

  row.names(tri)<-c("A","B","C")  #vertex labeling

  Avec= paste(A, collapse=","); Bvec= paste(B, collapse=","); Cvec= paste(C, collapse=",");
  typ<-paste("Closest Points in CM-Vertex Regions of the Triangle with Vertices A=(",Avec,"), B=(",Bvec,"), and C=(",Cvec,") \n to the Opposite Edges",sep="")
  description<-"Closest Points in CM-Vertex Regions of the Triangle to its Edges \n (Row i corresponds to vertex region i for i=1,2,3)"
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the Edges and the Closest Points to the Edges in CM-Vertex Regions \n (i-th entry corresponds to vertex region i for i=1,2,3)"
  main.txt<-"Closest Points in CM-Vertex Regions \n to the Opposite Edges"

  Dis<-c(ifelse(is.numeric(distA),min(distA),NA),ifelse(is.numeric(distB),min(distB),NA),
         ifelse(is.numeric(distC),min(distC),NA))
  #distances of the closest points to the edges in the respective vertex regions
  Regs<-list(vr1=rbind(A,D3,Cent,D2), #regions inside the triangles
             vr2=rbind(B,D1,Cent,D3),
             vr3=rbind(C,D2,Cent,D1)
  )

  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3")  #regions names
  supp.type = "Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=ce,  #closest points to edges in each associated vertex region
    X=Xp, num.points=n, #data points and its size
    ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (triangle here)
    cent=Cent, cent.name=Cname, #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points in vertex regions to the corresponding edges
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points in a data set to edges in each \eqn{CC}-vertex region in a triangle
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to edge \eqn{j} in \eqn{CC}-vertex region \eqn{j} for \eqn{j=1,2,3}
#' in the triangle, \code{tri}\eqn{=T(A,B,C)}, where \eqn{CC} stands for circumcenter. Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3},
#' and corresponding edge labels are \eqn{BC=1}, \eqn{AC=2}, and \eqn{AB=3}.
#' Function yields \code{NA} if there are no data points in a \eqn{CC}-vertex region.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010;textual}{pcds}).
#'
#' @param Xp A set of 2D points representing the set of data points.
#' @param tri A \eqn{3 \times 2} matrix with each row representing a vertex of the triangle.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, and \eqn{C=3} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances to Edges in the Respective
#' CC-Vertex Regions"}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to edges in the respective vertex region.}
#' \item{ind.ext}{Indices of the extrema points,\code{ext}.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \code{tri}}
#' \item{cent}{The center point used for construction of vertex regions}
#' \item{ncent}{Name of the center, \code{cent}, it is \code{"CC"} for this function}
#' \item{regions}{Vertex regions inside the triangle, \code{tri}, provided as a \code{list}}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1"}, \code{"vr=2"}, and \code{"vr=3"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{tri}}
#' \item{dist2ref}{Distances of closest points in the vertex regions to corresponding edges}
#'
#' @seealso \code{\link{cl2edgesVRbas.tri}}, \code{\link{cl2edgesCM_VR}}, \code{\link{cl2edgesM_VR}},
#' and \code{\link{cl2edges.std.tri}}
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
#' n<-20  #try also n<-100
#' set.seed(1)
#' Xp<-runif.tri(n,Tr)$g
#'
#' Ext<-cl2edgesCC_VR(Xp,Tr)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' cl2e<-Ext
#'
#' CC<-circ.cent.tri(Tr);
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#'
#' Xlim<-range(Tr[,1],Xp[,1],CC[1])
#' Ylim<-range(Tr[,2],Xp[,2],CC[2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(Tr,asp=1,pch=".",xlab="",ylab="",
#' main="Closest Points in CC-Vertex Regions \n to the Opposite Edges",
#' axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Tr)
#'
#' xc<-Tr[,1]+c(-.02,.02,.02)
#' yc<-Tr[,2]+c(.02,.02,.04)
#' txt.str<-c("A","B","C")
#' text(xc,yc,txt.str)
#'
#' points(Xp,pch=1,col=1)
#' L<-matrix(rep(CC,3),ncol=2,byrow=TRUE); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#' points(cl2e$ext,pch=3,col=2)
#'
#' txt<-rbind(CC,Ds)
#' xc<-txt[,1]+c(-.04,.04,-.03,0)
#' yc<-txt[,2]+c(-.05,.04,.06,-.08)
#' txt.str<-c("CC","D1","D2","D3")
#' text(xc,yc,txt.str)
#' }
#'
#' @export cl2edgesCC_VR
cl2edgesCC_VR <- function(Xp,tri)
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

  n<-nrow(Xp)

  A<-tri[1,]; B<-tri[2,]; C<-tri[3,]
  Cent<-circ.cent.tri(tri) ; Cname<-"CC"
  D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;

  clBC<-NA; clAC<-NA; clAB<-NA
  VRdt<-rel.verts.triCC(Xp,tri)$rv
  dtA<-matrix(Xp[VRdt==1 & !is.na(VRdt),],ncol=2); #this is data points in CC-vertex region A
  dtB<-matrix(Xp[VRdt==2 & !is.na(VRdt),],ncol=2);
  dtC<-matrix(Xp[VRdt==3 & !is.na(VRdt),],ncol=2)

  ind.dtA<-which(VRdt==1 & !is.na(VRdt)); #indices of the data points in CC-vertex region A
  ind.dtB<-which(VRdt==2 & !is.na(VRdt));
  ind.dtC<-which(VRdt==3 & !is.na(VRdt))

  distA<-distB<-distC<-vector()
  nA<-nrow(dtA); nB<-nrow(dtB); nC<-nrow(dtC);
  clBC.ind = clAC.ind = clAB.ind = NA

  if (nA>0)
  {
    for (i in 1:nA)
    {distA<-c(distA,dist.point2line(dtA[i,],B,C)$dis)};
    clBC<-dtA[distA==min(distA),]  #closest to edge BC in vertex region 1
    clBC.ind = ind.dtA[distA==min(distA)] #index of the closest point
  }

  if (nB>0)
  {
    for (i in 1:nB)
    {distB<-c(distB,dist.point2line(dtB[i,],A,C)$dis)};
    clAC<-dtB[distB==min(distB),]
    clAC.ind = ind.dtB[distB==min(distB)]
  }

  if (nC>0)
  {
    for (i in 1:nC)
    {distC<-c(distC,dist.point2line(dtC[i,],A,B)$dis)};
    clAB<-dtC[distC==min(distC),]
    clAB.ind = ind.dtC[distC==min(distC)]
  }

  ce<-rbind(clBC,clAC,clAB) #closest points to opposite edges in the vertex regions
  row.names(ce) = NULL
  ce.ind = c(clBC.ind,clAC.ind,clAB.ind) #their indices

  Dis<-c(ifelse(is.numeric(distA),min(distA),NA),ifelse(is.numeric(distB),min(distB),NA),
             ifelse(is.numeric(distC),min(distC),NA))

  row.names(tri)<-c("A","B","C")  #vertex labeling
  txt<-"Edge labels are AB=3, BC=1, and AC=2 (correspond to row number in ce)"

  row.names(tri)<-c("A","B","C")  #vertex labeling

  Avec= paste(A, collapse=","); Bvec= paste(B, collapse=","); Cvec= paste(C, collapse=",");
  typ<-paste("Closest Points in CC-Vertex Regions of the Triangle with Vertices A=(",Avec,"), B=(",Bvec,"), and C=(",Cvec,") \n to the Opoosite Edges",sep="")
  description<-"Closest Points in CC-Vertex Regions of the Triangle to its Edges \n (Row i corresponds to vertex region i for i=1,2,3)"
  txt1<-"Vertex labels are A=1, B=2, and C=3 (correspond to row number in Extremum Points)"
  txt2<-"Distances between the Edges and the Closest Points to the Edges in CC-Vertex Regions \n (i-th entry corresponds to vertex region i for i=1,2,3)"
  main.txt<-"Closest Points in CC-Vertex Regions \n to the Opposite Edges"

  #distances of the closest points to the edges in the respective vertex regions
  if (in.triangle(Cent,tri,boundary = FALSE)$i)
  {
    Regs<-list(vr1=rbind(A,D3,Cent,D2), #regions inside the triangles
               vr2=rbind(B,D1,Cent,D3),
               vr3=rbind(C,D2,Cent,D1))
    Reg.names<-c("vr=1","vr=2","vr=3")  #regions names
  } else
  { a1<-A[1]; a2<-A[2]; b1<-B[1]; b2<-B[2]; c1<-C[1]; c2<-C[2];
  dAB<-Dist(A,B); dBC<-Dist(B,C); dAC<-Dist(A,C); max.dis<-max(dAB,dBC,dAC)
  if (dAB==max.dis)
  {
    L1<-c((1/2)*(a1^3-a1^2*b1+a1*a2^2-2*a1*a2*b2+2*a1*b2*c2-a1*c1^2-a1*c2^2+a2^2*b1-2*a2*b1*c2+b1*c1^2+b1*c2^2)/(a1^2-a1*b1-a1*c1+a2^2-a2*b2-a2*c2+b1*c1+b2*c2), (1/2)*(a1^2*a2+a1^2*b2-2*a1*a2*b1-2*a1*b2*c1+a2^3-a2^2*b2+2*a2*b1*c1-a2*c1^2-a2*c2^2+b2*c1^2+b2*c2^2)/(a1^2-a1*b1-a1*c1+a2^2-a2*b2-a2*c2+b1*c1+b2*c2))
    L2<-c((1/2)*(a1*b1^2-a1*b2^2+2*a1*b2*c2-a1*c1^2-a1*c2^2+2*a2*b1*b2-2*a2*b1*c2-b1^3-b1*b2^2+b1*c1^2+b1*c2^2)/(a1*b1-a1*c1+a2*b2-a2*c2-b1^2+b1*c1-b2^2+b2*c2), (1/2)*(2*a1*b1*b2-2*a1*b2*c1-a2*b1^2+2*a2*b1*c1+a2*b2^2-a2*c1^2-a2*c2^2-b1^2*b2-b2^3+b2*c1^2+b2*c2^2)/(a1*b1-a1*c1+a2*b2-a2*c2-b1^2+b1*c1-b2^2+b2*c2))
    Regs<-list(vr1=rbind(A,L1,D2), #regions inside the triangles
               vr2=rbind(B,D1,L2),
               vr3=rbind(C,D2,L1,L2,D1),
               r4=rbind(Cent,L1,L2))  #only r4 is outside the triangle
  }

  if (dBC==max.dis)
  {
    L1<-c((1/2)*(a1^2*b1-a1^2*c1+a2^2*b1-a2^2*c1-2*a2*b1*c2+2*a2*b2*c1-b1^3+b1^2*c1-b1*b2^2+2*b1*b2*c2-b2^2*c1)/(a1*b1-a1*c1+a2*b2-a2*c2-b1^2+b1*c1-b2^2+b2*c2), (1/2)*(a1^2*b2-a1^2*c2+2*a1*b1*c2-2*a1*b2*c1+a2^2*b2-a2^2*c2-b1^2*b2-b1^2*c2+2*b1*b2*c1-b2^3+b2^2*c2)/(a1*b1-a1*c1+a2*b2-a2*c2-b1^2+b1*c1-b2^2+b2*c2))
    L2<-c((1/2)*(a1^2*b1-a1^2*c1+a2^2*b1-a2^2*c1-2*a2*b1*c2+2*a2*b2*c1-b1*c1^2+b1*c2^2-2*b2*c1*c2+c1^3+c1*c2^2)/(a1*b1-a1*c1+a2*b2-a2*c2-b1*c1-b2*c2+c1^2+c2^2), (1/2)*(a1^2*b2-a1^2*c2+2*a1*b1*c2-2*a1*b2*c1+a2^2*b2-a2^2*c2-2*b1*c1*c2+b2*c1^2-b2*c2^2+c1^2*c2+c2^3)/(a1*b1-a1*c1+a2*b2-a2*c2-b1*c1-b2*c2+c1^2+c2^2))
    Regs<-list(vr1=rbind(A,D3,L1,L2,D2), #regions inside the triangles
               vr2=rbind(B,L1,D3),
               vr3=rbind(C,D2,L2),
               r4=rbind(Cent,L1,L2))  #only r4 is outside the triangle
  }

  if (dAC==max.dis)
  {
    L1<-c((1/2)*(a1^3-a1^2*c1+a1*a2^2-2*a1*a2*c2-a1*b1^2-a1*b2^2+2*a1*b2*c2+a2^2*c1-2*a2*b2*c1+b1^2*c1+b2^2*c1)/(a1^2-a1*b1-a1*c1+a2^2-a2*b2-a2*c2+b1*c1+b2*c2), (1/2)*(a1^2*a2+a1^2*c2-2*a1*a2*c1-2*a1*b1*c2+a2^3-a2^2*c2-a2*b1^2+2*a2*b1*c1-a2*b2^2+b1^2*c2+b2^2*c2)/(a1^2-a1*b1-a1*c1+a2^2-a2*b2-a2*c2+b1*c1+b2*c2))
    L2<-c((1/2)*(a1*b1^2+a1*b2^2-2*a1*b2*c2-a1*c1^2+a1*c2^2+2*a2*b2*c1-2*a2*c1*c2-b1^2*c1-b2^2*c1+c1^3+c1*c2^2)/(a1*b1-a1*c1+a2*b2-a2*c2-b1*c1-b2*c2+c1^2+c2^2), (1/2)*(2*a1*b1*c2-2*a1*c1*c2+a2*b1^2-2*a2*b1*c1+a2*b2^2+a2*c1^2-a2*c2^2-b1^2*c2-b2^2*c2+c1^2*c2+c2^3)/(a1*b1-a1*c1+a2*b2-a2*c2-b1*c1-b2*c2+c1^2+c2^2))
    Regs<-list(vr1=rbind(A,D3,L1), #regions inside the triangles
               vr2=rbind(B,D1,L2,L1,D3),
               vr3=rbind(C,L2,D1),
               r4=rbind(Cent,L1,L2))  #only r4 is outside the triangle
  }
  Reg.names<-c("vr=1","vr=2","vr=3",NA)  #regions names
  }

  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  supp.type = "Triangle" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=ce, ind.ext= ce.ind, #closest points to edges in each associated vertex region, and their indices
    X=Xp, num.points=n, #data points and its size
    ROI=tri, supp.type = supp.type, # region of interest for X points, and its type (triangle here)
    cent=Cent, cent.name=Cname, #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points in vertex regions to the corresponding edges
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The closest points among a data set in the vertex regions to the respective faces in a tetrahedron
#'
#' @description
#' An object of class \code{"Extrema"}.
#' Returns the closest data points among the data set, \code{Xp}, to face \eqn{i} in \code{M}-vertex region \eqn{i} for \eqn{i=1,2,3,4}
#' in the tetrahedron \eqn{th=T(A,B,C,D)}. Vertex labels are \eqn{A=1}, \eqn{B=2}, \eqn{C=3}, and \eqn{D=4} and corresponding face
#' labels are \eqn{BCD=1}, \eqn{ACD=2}, \eqn{ABD=3}, and \eqn{ABC=4}.
#'
#' Vertex regions are based on center \code{M} which can be the center
#' of mass (\code{"CM"}) or circumcenter (\code{"CC"}) of \code{th}.
#'
#' @param Xp A set of 3D points representing the set of data points.
#' @param th A \eqn{4 \times 3} matrix with each row representing a vertex of the tetrahedron.
#' @param M The center to be used in the construction of the vertex regions in the tetrahedron, \code{th}.
#' Currently it only takes \code{"CC"} for circumcenter and \code{"CM"} for center of mass; default=\code{"CM"}.
#'
#' @return A \code{list} with the elements
#' \item{txt1}{Vertex labels are \eqn{A=1}, \eqn{B=2}, \eqn{C=3}, and \eqn{D=4} (correspond to row number in Extremum Points).}
#' \item{txt2}{A short description of the distances as \code{"Distances from Closest Points to Faces ..."}.}
#' \item{type}{Type of the extrema points}
#' \item{desc}{A short description of the extrema points}
#' \item{mtitle}{The \code{"main"} title for the plot of the extrema}
#' \item{ext}{The extrema points, here, closest points to faces in the respective vertex region.}
#' \item{ind.ext}{The data indices of extrema points, \code{ext}.}
#' \item{X}{The input data, \code{Xp}, can be a \code{matrix} or \code{data frame}}
#' \item{num.points}{The number of data points, i.e., size of \code{Xp}}
#' \item{supp}{Support of the data points, here, it is \code{th}}
#' \item{cent}{The center point used for construction of vertex regions, it is circumcenter of center of mass
#' for this function}
#' \item{ncent}{Name of the center, it is circumcenter \code{"CC"} or  center of mass \code{"CM"} for this function.}
#' \item{regions}{Vertex regions inside the tetrahedron \code{th} provided as a list.}
#' \item{region.names}{Names of the vertex regions as \code{"vr=1","vr=2","vr=3","vr=4"}}
#' \item{region.centers}{Centers of mass of the vertex regions inside \code{th}.}
#' \item{dist2ref}{Distances from closest points in each vertex region to the corresponding face.}
#'
#' @seealso \code{\link{fr2vertsCC_VR}}, \code{\link{fr2edgesCM_ERstd.tri}}, \code{\link{kfr2vertsCC_VRbas.tri}} and \code{\link{kfr2vertsCC_VR}}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' A<-c(0,0,0); B<-c(1,0,0); C<-c(1/2,sqrt(3)/2,0); D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
#' set.seed(1)
#' tetra<-rbind(A,B,C,D)+matrix(runif(12,-.25,.25),ncol=3)
#' n<-10  #try also n<-20
#' Cent<-"CC"  #try also "CM"
#'
#' n<-20  #try also n<-100
#' Xp<-runif.tetra(n,tetra)$g  #try also Xp<-cbind(runif(n),runif(n),runif(n))
#'
#' Ext<-cl2facesVRtetra(Xp,tetra,Cent)
#' Ext
#' summary(Ext)
#' plot(Ext)
#'
#' clf<-Ext$ext
#'
#' if (Cent=="CC") {M<-circ.cent.tetra(tetra)}
#' if (Cent=="CM") {M<-apply(tetra,2,mean)}
#'
#' Xlim<-range(tetra[,1],Xp[,1],M[1])
#' Ylim<-range(tetra[,2],Xp[,2],M[2])
#' Zlim<-range(tetra[,3],Xp[,3],M[3])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#' zd<-Zlim[2]-Zlim[1]
#'
#' plot3D::scatter3D(Xp[,1],Xp[,2],Xp[,3], phi =0,theta=40, bty = "g",
#' main="Closest Pointsin CC-Vertex Regions \n to the Opposite Faces",
#' xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05), zlim=Zlim+zd*c(-.05,.05),
#'           pch = 20, cex = 1, ticktype = "detailed")
#' #add the vertices of the tetrahedron
#' plot3D::points3D(tetra[,1],tetra[,2],tetra[,3], add=TRUE)
#' L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
#' plot3D::points3D(clf[,1],clf[,2],clf[,3], pch=4,col="red", add=TRUE)
#'
#' plot3D::text3D(tetra[,1],tetra[,2],tetra[,3], labels=c("A","B","C","D"), add=TRUE)
#'
#' #for center of mass use #Cent<-apply(tetra,2,mean)
#' D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
#' L<-rbind(D1,D2,D3,D4,D5,D6); R<-rbind(M,M,M,M,M,M)
#' plot3D::segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
#' }
#'
#' @export cl2facesVRtetra
cl2facesVRtetra <- function(Xp,th,M="CM")
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (is.point(Xp,3))
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

  if (length(M) > 1 || sum(M==c("CM","CC"))==0)
    stop("M must be one of \"CC\", \"CM\"")

  y1<-A<-th[1,]; y2<-B<-th[2,]; y3<-C<-th[3,]; y4<-D<-th[4,];
  ifelse(identical(M,"CC"), {Cent<-circ.cent.tetra(th); cent.name<-"CC"},
         {Cent<-apply(th, 2, mean); cent.name<-"CM"})
  D1<-(y1+y2)/2; D2<-(y1+y3)/2; D3<-(y1+y4)/2; D4<-(y2+y3)/2; D5<-(y2+y4)/2; D6<-(y3+y4)/2;
  Ds<-rbind(D1,D2,D3,D4,D5,D6)

  mdt<-c(dist.point2plane(y1,y2,y3,y4)$d,dist.point2plane(y2,y1,y3,y4)$d,dist.point2plane(y3,y1,y2,y4)$d,dist.point2plane(y4,y1,y2,y3)$d);
  #distances from each vertex to the opposite face

  U<-matrix(NA,nrow=4,ncol=3); # closest data points in vertex regions to their respective (opposite to the vertex) faces
  ind.ext<-rep(NA,4)  # data indices of the closest points in vertex regions to their respective (opposite to the vertex) faces

  Xp<-matrix(Xp,ncol=3)
  n<-nrow(Xp);

  for (i in 1:n)
  {ifelse(identical(M,"CC"),rv<-rv.tetraCC(Xp[i,],th)$rv,rv<-rv.tetraCM(Xp[i,],th)$rv);
    if (!is.na(rv))
    {
      if (rv==1)
      {d1<-dist.point2plane(Xp[i,],y2,y3,y4)$d;
      if (d1<=mdt[1]) {mdt[1]<-d1; U[1,]<-Xp[i,]; ind.ext[1]=i};
      } else if (rv==2)
      {d2<-dist.point2plane(Xp[i,],y1,y3,y4)$d;
      if (d2<=mdt[2]) {mdt[2]<-d2; U[2,]<-Xp[i,]; ind.ext[2]=i}
      } else if (rv==3)
      {d3<-dist.point2plane(Xp[i,],y1,y2,y4)$d;
      if (d3<=mdt[3]) {mdt[3]<-d3; U[3,]<-Xp[i,]; ind.ext[3]=i}
      } else {
        d4<-dist.point2plane(Xp[i,],y1,y2,y3)$d;
        if (d4<=mdt[4]) {mdt[4]<-d4; U[4,]<-Xp[i,]; ind.ext[4]=i}
      }
    }
  }

  row.names(th)<-c("A","B","C","D")  #vertex labeling

  Avec= paste(round(y1,2), collapse=","); Bvec= paste(round(y2,2), collapse=",");
  Cvec= paste(round(y3,2), collapse=","); Dvec= paste(round(y4,2), collapse=",");
  typ<-paste("Closest Points in ",cent.name,"-Vertex Regions of the Tetrahedron with Vertices A=(",Avec,"), B=(",Bvec,"), C=(",Cvec,"), and D=(",Dvec,") to the Opposite Faces",sep="")
  description<-paste("Closest Points in ",cent.name,"-Vertex Regions of the Tetrahedron to its Faces\n (Row i corresponds to face i for i=1,2,3,4)",sep="")
  txt1<-"Vertex labels are A=1, B=2, C=3, and D=4 (correspond to row number in Extremum Points)"
  txt2<-paste("Distances between Faces and the Closest Points to the Faces in ",cent.name,"-Vertex Regions\n (i-th entry corresponds to vertex region i for i=1,2,3,4)",sep="")
  main.txt<-paste("Closest Points in ",cent.name,"-Vertex Regions \n of the Tetrahedron to its Faces",sep="")

  Dis<-c(ifelse(!is.na(U[1,1]),mdt[1],NA),ifelse(!is.na(U[2,1]),mdt[2],NA),
             ifelse(!is.na(U[3,1]),mdt[3],NA),ifelse(!is.na(U[4,1]),mdt[4],NA))
  #distances of the closest points to the faces in the respective vertex regions
  Regs<-list(vr1=rbind(A,D1,D2,D3,Cent), #regions inside the triangles
             vr2=rbind(B,D1,D4,D5,Cent),
             vr3=rbind(C,D4,D2,D6,Cent),
             vr4=rbind(D,D3,D5,D6,Cent)
  )
  Reg.Cent<-vector()
  for (i in 1:length(Regs))
  { Reg.Cent<-rbind(Reg.Cent,apply(Regs[[i]],2,mean))}

  Reg.names<-c("vr=1","vr=2","vr=3","vr=4")  #regions names

  supp.type = "Tetrahedron" #name of the support

  res<-list(
    txt1=txt1, txt2=txt2,
    type=typ, desc=description,
    mtitle=main.txt, #main label in the plot
    ext=U, ind.ext=ind.ext, #ext: closest points to faces in each associated vertex region, and ind.ext is their data indices
    X=Xp, num.points=n, #data points and its size
    ROI=th, supp.type = supp.type, # region of interest for X points, and its type (tetrahedron here)
    cent=Cent, cent.name=cent.name, #center and center name
    regions=Regs, region.names=Reg.names, region.centers=Reg.Cent,
    dist2ref=Dis #distances of closest points to edges in the respective vertex region
  )

  class(res)<-"Extrema"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

# funsRankOrderTe
#'
#' @title Returns the ranks and orders of points in decreasing distance to the edges of the triangle
#'
#' @description
#' Two functions: \code{rank.d2e.Te} and \code{order.d2e.Te}.
#'
#' \code{rank.d2e.Te} finds the ranks of the distances of points in data, \code{Xp}, to the edges of the standard equilateral triangle
#' \eqn{T_e=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#'
#' \code{dec} is a logical argument, default is \code{TRUE}, so the ranks are for decreasing distances, if \code{FALSE} it will be
#' in increasing distances.
#'
#' \code{order.d2e.Te} finds the orders of the distances of points in data, \code{Xp}, to the edges of \eqn{T_e}. The arguments are
#' as in \code{rank.d2e.Te}.
#'
#' @param Xp A set of 2D points representing the data set in which ranking in terms of the distance
#' to the edges of \eqn{T_e} is performed.
#' @param dec A logical argument indicating the how the ranking will be performed. If \code{TRUE},
#' the ranks are for decreasing distances,and if \code{FALSE} they will be in increasing distances,
#' default is \code{TRUE}.
#'
#' @return A \code{list} with two elements
#' \item{distances}{Distances from data points to the edges of \eqn{T_e}}
#' \item{dist.rank}{The ranks of the data points in decreasing distances to the edges of \eqn{T_e}}
#'
#' @name funsRankOrderTe
NULL
#'
#' @rdname funsRankOrderTe
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \dontrun{
#' #Examples for rank.d2e.Te
#' n<-10
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points
#'
#' dec.dist<-rank.d2e.Te(Xp)
#' dec.dist
#' dec.dist.rank<-dec.dist[[2]] #the rank of distances to the edges in decreasing order
#' dec.dist.rank
#'
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.0,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp,pch=".")
#' text(Xp,labels = factor(dec.dist.rank) )
#'
#' inc.dist<-rank.d2e.Te(Xp,dec = FALSE)
#' inc.dist
#' inc.dist.rank<-inc.dist[[2]] #the rank of distances to the edges in increasing order
#' inc.dist.rank
#' dist<-inc.dist[[1]] #distances to the edges of the std eq. triangle
#' dist
#'
#' plot(A,pch=".",xlab="",ylab="",xlim=Xlim,ylim=Ylim)
#' polygon(Te)
#' points(Xp,pch=".",xlab="",ylab="", main="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' text(Xp,labels = factor(inc.dist.rank))
#' }
#'
#' @export
rank.d2e.Te <- function(Xp,dec=TRUE)
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

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)

  n<-nrow(Xp)
  dist.edge<-rep(NA,n)
  for (i in 1:n)
  { pnt<-as.vector(Xp[i,])
  if (in.triangle(pnt,Te,boundary = TRUE)$in.tri)
  {
    x<-pnt[1]; y<-pnt[2];
    dist.edge[i]<-min(y,-0.5*y+0.866025404*x,-0.5*y+0.8660254040-0.866025404*x)
  }
  }

  ifelse(dec==TRUE,ranks<-rank(-dist.edge),ranks<-rank(dist.edge))
  ranks[is.na(dist.edge)]<-NA

  list(distances=dist.edge,
       dist.rank=ranks)
} #end of the function
#'
#' @rdname funsRankOrderTe
#'
#' @examples
#' \dontrun{
#' #Examples for order.d2e.Te
#' n<-10
#' set.seed(1)
#' Xp<-runif.std.tri(n)$gen.points  #try also Xp<-cbind(runif(n),runif(n))
#'
#' dec.dist<-order.d2e.Te(Xp)
#' dec.dist
#' dec.dist.order<-dec.dist[[2]] #the order of distances to the edges in decreasing order
#' dec.dist.order
#'
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C);
#'
#' Xlim<-range(Te[,1])
#' Ylim<-range(Te[,2])
#' xd<-Xlim[2]-Xlim[1]
#' yd<-Ylim[2]-Ylim[1]
#'
#' plot(A,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.01,.01),ylim=Ylim+yd*c(-.01,.01))
#' polygon(Te)
#' points(Xp,pch=".")
#' text(Xp[dec.dist.order,],labels = factor(1:n) )
#'
#' inc.dist<-order.d2e.Te(Xp,dec = FALSE)
#' inc.dist
#' inc.dist.order<-inc.dist[[2]] #the order of distances to the edges in increasing order
#' inc.dist.order
#' dist<-inc.dist[[1]] #distances to the edges of the std eq. triangle
#' dist
#' dist[inc.dist.order] #distances in increasing order
#'
#' plot(A,pch=".",xlab="",ylab="",xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
#' polygon(Te)
#' points(Xp,pch=".")
#' text(Xp[inc.dist.order,],labels = factor(1:n))
#' }
#'
#' @export
order.d2e.Te <- function(Xp,dec=TRUE)
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

  A<-c(0,0); B<-c(1,0); C<-c(0.5,sqrt(3)/2); Te<-rbind(A,B,C)

  n<-nrow(Xp)
  dist.edge<-rep(NA,n)
  for (i in 1:n)
  {  pnt<-as.vector(Xp[i,])
  if (in.triangle(pnt,Te,boundary = TRUE)$in.tri)
  {
    x<-pnt[1]; y<-pnt[2];
    dist.edge[i]<-min(y,-0.5*y+0.866025404*x,-0.5*y+0.8660254040-0.866025404*x)
  }
  }
  ifelse(dec==TRUE,orders<-order(dist.edge,decreasing=dec),orders<-order(dist.edge))
  nint<-sum(!is.na(dist.edge))
  orders[-(1:nint)]<-NA

  list(distances=dist.edge,dist.order=orders)
} #end of the function



