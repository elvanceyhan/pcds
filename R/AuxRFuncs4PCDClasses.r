#AuxRFuncs4PCDClasses.r
###############################################
#Auxiliary functions for class Lines
###############################################
print.Lines <- function(x, ...)
{
  if (!inherits(x, "Lines"))
    stop("x must be of class \"Lines\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients of the line (for the form: y = slope * x + intercept ) \n")
  print(c(x$slope,x$intercept))
}

########################
summary.Lines <- function(object, ...)
{
  if (!inherits(object, "Lines"))
    stop("object must be of class \"Lines\"")
  
  eqn <- object$equation
  names(eqn)<-c()
  pnts<-object$points
  xv<-object$x
  yv<-object$y
  nx<-min(6,length(xv))
  ny<-min(6,length(yv))
  
  res <- list(txt=object$txt,
              call=object$call,
              points=pnts,
              xvec=xv[1:nx],
              yvec=yv[1:ny],
              coefficients=c(object$slope,object$intercept),
              eqn=eqn)
  class(res) <- "summary.Lines"
  res
}

########################

print.summary.Lines <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nDefining Points\n")
  print(x$points)
  cat("\n Selected x points (first row) and estimated y points (second row) that fall on the Line
      (first 6 or fewer are printed at each row) \n")
  print(x$xvec)
  print(x$yvec)
  
  cat("\nEquation of the", x$txt,"\n")
  print(x$eqn)
  
  cat("\nCoefficients of the line (when the line is in the form: y = slope * x + intercept ) \n")
  print(x$coefficients)
}

########################
plot.Lines<-function (x, asp=NA, xlab="x",ylab="y",...)
{
  if (!inherits(x, "Lines"))
    stop("x must be of class \"Lines\"")
  
  A<-x$points[1,]
  B<-x$points[2,]
  pts<-x$points
  vert<-x$vert
  eqn.lab<-x$eqnlabel
  
  xl = range(x$x,pts)
  
  xf<-(xl[2]-xl[1])*.1
  xp<-seq(xl[1]-xf,xl[2]+xf,l=100)
  
  Sl<-x$slope
  if (abs(Sl)==Inf)
  {yp<-NA
  Ylim<-range(0,pts[,2])
  } else
  {
    yp<-x$intercept+Sl*xp
    Ylim<-range(yp,pts[,2])
  }
  
  Xlim<-range(xp,pts[,1])
  
  xd<-Xlim[2]-Xlim[1]
  yd<-Ylim[2]-Ylim[1]
  pf<-c(xd,-yd)*.025
  ptf<-cbind(rep(pf[1],nrow(pts)),rep(pf[2],nrow(pts)))
  if (abs(Sl)==Inf)
  {
    plot(xp,rep(0,length(xp)),type="n",xlab=xlab,ylab=ylab,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),main=x$mtxt)
    abline(v=vert) 
    abline(h=0,lty=2) 
    text(rbind(c(vert,0))+pf*.01, col=1,eqn.lab)
  } else
  {
    plot(xp,yp,asp=asp,type="l",xlab=xlab,ylab=ylab,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),main=x$mtxt)
    text(rbind(c(median(xp),median(yp)))+pf*3, col=1,eqn.lab)
  }
  points(pts)
  text(pts+ptf,row.names(pts))
  segments(A[1],A[2], B[1], B[2], lty=1,col=2,lwd=2)
}

###############################################
#Auxiliary functions for class Lines3D
###############################################
print.Lines3D <- function(x, ...)
{
  if (!inherits(x, "Lines3D"))
    stop("x must be of class \"Lines3D\"")
  
  cat("Call:\n")
  print(x$call)
  vnames<-x$vec.names
  cat("\nCoefficients of the parameterized line passing through initial point", vnames[1],"=(x0,y0,z0) in the direction of",vnames[2],"=(a,b,c) (for the form: x=x0 + a*t, y=y0 + b*t, and z=z0 + c*t ) \n")
  print(x$vecs)
}

########################
summary.Lines3D <- function(object, ...)
{
  if (!inherits(object, "Lines3D"))
    stop("object must be of class \"Lines3D\"")
  
  eqn <- object$equation
  names(eqn)<-c()
  vecs<-vecs2<-object$vecs
  vnames<-object$vec.names
  row.names(vecs2)<-vnames
  xv<-object$x
  yv<-object$y
  zv<-object$z
  nx<-min(6,length(xv))
  ny<-min(6,length(yv))
  nz<-min(6,length(zv))
  
  res <- list(txt=object$txt,
              call=object$call,
              vectors=vecs,
              xvec=xv[1:nx],
              yvec=yv[1:ny],
              zvec=zv[1:nz],
              coefficients=vecs2,
              eqn=eqn)
  class(res) <- "summary.Lines3D"
  res
}

########################

print.summary.Lines3D <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nDefining Vectors\n")
  print(x$vectors)
  cat("\n Estimated x points (first row), y points (second row), and z points (third row)
   that fall on the Line
      (first 6 or fewer are printed at each row) \n")
  print(x$xvec)
  print(x$yvec)
  print(x$zvec)
  
  cat("\nEquation of", x$txt,"\n")
  print(x$eqn)
  
  Vnames<-row.names(x$coefficients)
  cat("\nCoefficients of the parameterized line passing through initial point", Vnames[1],"=(x0,y0,z0) in the direction of ",Vnames[2],"=(a,b,c) (for the form: x=x0 + a*t, y=y0 + b*t, and z=z0 + c*t ) \n")
  print(x$coefficients)
}

########################
plot.Lines3D<-function (x,asp=NA,xlab="x",ylab="y",zlab="z",...)
{
  if (!inherits(x, "Lines3D"))
    stop("x must be of class \"Lines3D\"")
  
  A<-x$vecs[1,]
  Bd<-x$vecs[2,]
  Pts<-x$pts
  
  xc<-x$x
  yc<-x$y
  zc<-x$z
  
  tr<-range(x$tsq); 
  tf<-(tr[2]-tr[1])*.1
  
  zr<-range(zc)
  zf<-(zr[2]-zr[1])*.2
  Bv<-Bd*tf*5
  
  Xlim<-range(xc,Pts[,1])
  Ylim<-range(yc,Pts[,2])
  Zlim<-range(zc,Pts[,3])
  
  xd<-Xlim[2]-Xlim[1]
  yd<-Ylim[2]-Ylim[1]
  zd<-Zlim[2]-Zlim[1]
  
  Dr<-A+min(x$tsq)*Bd
  
  if (x$vec.names[2]=="normal vector")
  {
   A0<-Pts[2,]; B0<-Pts[3,]; C0<-Pts[4,];
   pts<-rbind(A0,B0,C0)
   xr<-range(pts[,1],xc); yr<-range(pts[,2],yc)
   xf<-(xr[2]-xr[1])*.1 #how far to go at the lower and upper ends in the x-coordinate
   yf<-(yr[2]-yr[1])*.1 #how far to go at the lower and upper ends in the y-coordinate
   xs<-seq(xr[1]-xf,xr[2]+xf,l=100)
   ys<-seq(yr[1]-yf,yr[2]+yf,l=100)
    
   plABC<-Plane(A0,B0,C0,xs,ys)
   z.grid<-plABC$z
    
   Zlim<-range(zc,z.grid,Pts[,3])
   zd<-Zlim[2]-Zlim[1]
    
   Bv<- -Bd*zf*5
   
   Dr<-(A0+B0+C0)/3
   
   persp3D(z = z.grid, x = xs, y = ys, asp=1,theta =225, phi = 30, col="lightblue", ticktype = "detailed",xlab=xlab,ylab=ylab,zlab=zlab,
            xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.05,.05),main=x$mtxt) #plane spanned by points A, B, C
   lines3D(xc, yc, zc, bty = "g",pch = 20, cex = 2, ticktype = "detailed",add=T)
   if (!is.null(Pts))
   {points3D(Pts[,1],Pts[,2],Pts[,3],add=T)
     text3D(Pts[,1],Pts[,2],Pts[,3],labels=x$pnames,add=T)}
   arrows3D(A[1],A[2],A[3]-zf,A[1],A[2],A[3],lty=2, add=T)
   text3D(A[1],A[2],A[3]-zf,labels="initial point",add=T)
   text3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+zf+Bv[3]/2,labels=x$vec.names[2],add=T)
   arrows3D(Dr[1],Dr[2],Dr[3]+zf,Dr[1]+Bv[1],Dr[2]+Bv[2],Dr[3]+zf+Bv[3], add=T)
   text3D(A[1],A[2],A[3]+zf/2,labels=x$vec.names[1],add=T)
   arrows3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+3*zf+Bv[3]/2,Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+zf+Bv[3]/2,lty=2, add=T)
   
  } else
  {
  lines3D(xc, yc, zc, phi = 0, bty = "g",main=x$mtxt,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.1,.1),
           cex = 2, ticktype = "detailed")
    if (!is.null(Pts))
    {points3D(Pts[,1],Pts[,2],Pts[,3],add=T)
    text3D(Pts[,1],Pts[,2],Pts[,3],labels=x$pnames,add=T)}
    arrows3D(A[1],A[2],A[3]-zf,A[1],A[2],A[3],lty=2, add=T)
    text3D(A[1],A[2],A[3]-zf,labels="initial point",add=T)
    text3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+3*zf+Bv[3]/2,labels="direction vector",add=T)
    arrows3D(Dr[1],Dr[2],Dr[3]+zf,Dr[1]+Bv[1],Dr[2]+Bv[2],Dr[3]+zf+Bv[3], add=T)
    text3D(A[1],A[2],A[3]+zf/2,labels=x$vec.names[1],add=T)
    arrows3D(Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+3*zf+Bv[3]/2,Dr[1]+Bv[1]/2,Dr[2]+Bv[2]/2,Dr[3]+zf+Bv[3]/2,lty=2, add=T)
    }
}

###############################################
#Auxiliary functions for class Planes
###############################################
print.Planes <- function(x, ...)
{
  if (!inherits(x, "Planes"))
    stop("x must be of class \"Planes\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients of the Plane (in the z = A*x + B*y + C form):\n")
  print(c(x$coeff))
}

########################
summary.Planes <- function(object, ...)
{
  if (!inherits(object, "Planes"))
    stop("object must be of class \"Planes\"")
  
  eqn <- object$equation
  names(eqn)<-c()
  pnts<-object$points
  xv<-object$x
  yv<-object$y
  zv<-object$z
  
  nv<-min(6,length(xv),length(yv),length(zv))
  
  res <- list(txt=object$txt,
              call=object$call,
              points=pnts,
              xvec=xv[1:nv],
              yvec=yv[1:nv],
              zvec=zv[1:nv],
              coefficients=object$coeff,
              eqn=eqn)
  class(res) <- "summary.Planes"
  res
}

########################

print.summary.Planes <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nDefining Points\n")
  print(x$points)
  cat("\n Selected x and y points and estimated z points --- presented row-wise, respectively --- that fall on the Plane
      (first 6 or fewer are printed on each row) \n")
  print(x$xvec)
  print(x$yvec)
  print(x$zvec)
  
  cat("\nEquation of the", x$txt,"\n")
  print(x$eqn)
  
  cat("\nCoefficients of the Plane (in the z = A*x + B*y + C form):\n")
  print(x$coefficients)
}

########################
plot.Planes<-function (x, xlab="x",ylab="y",zlab="z",...)
{
  if (!inherits(x, "Planes"))
    stop("x must be of class \"Planes\"")
  
  A<-x$points[1,]
  B<-x$points[2,]
  C<-x$points[3,]
  pts<-x$points
  Cf<-x$coeff
  
  xl = range(x$x,pts[,1])
  yl = range(x$y,pts[,2])
  
  xf<-(xl[2]-xl[1])*.1
  yf<-(yl[2]-yl[1])*.1
  xp<-seq(xl[1]-xf,xl[2]+xf,l=100)
  yp<-seq(yl[1]-yf,yl[2]+yf,l=100)
  
  zgrid<-matrix(0,nrow=100,ncol=100)
  for (i in 1:100)
    for (j in 1:100)
    {
      zgrid[i,j]<-Cf[1]*xp[i]+Cf[2]*yp[j]+Cf[3]
    }
  
  xlim<-range(xp); ylim<-range(yp); zlim<-range(zgrid)
  
  xd<-xlim[2]-xlim[1]
  yd<-ylim[2]-ylim[1]
  zd<-zlim[2]-zlim[1]
  
  pxf<-xd*.025; pyf<-yd*.025; pzf<-zd*.025
  ptf<-cbind(rep(-pxf,nrow(pts)),rep(-pyf,nrow(pts)),rep(pzf,nrow(pts)))
  
  Pts<-pts+ptf
  Mn.pts<-c(mean(xp),mean(yp),mean(zgrid))
  
  library(plot3D)
  persp3D(z = zgrid, x = xp, y = yp, theta =225, phi = 30, ticktype = "detailed",xlab=xlab,ylab=ylab,zlab=zlab,
          xlim=xlim+xd*c(-.05,.05),ylim=ylim+yd*c(-.05,.05),zlim=zlim+zd*c(-.05,.05),main=x$txt) #plane spanned by points A, B, C
  points3D(Pts[,1],Pts[,2],Pts[,3], add = TRUE) #add the defining points
  text3D(Pts[,1],Pts[,2],Pts[,3], row.names(pts),add=TRUE)
  text3D(Mn.pts[1],Mn.pts[2],Mn.pts[3]+zd*.35,x$equation2,add=TRUE)
  polygon3D(Pts[1:3,1],Pts[1:3,2],Pts[1:3,3],add=TRUE)
  
}

###############################################
#Auxiliary functions for class TriLines
###############################################
print.TriLines <- function(x, ...)
{
  if (!inherits(x, "TriLines"))
    stop("x must be of class \"TriLines\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nEquation:\n")
  Eqn<-x$equation; names(Eqn)<-c()
  print(Eqn)
}


########################
summary.TriLines <- function(object, ...)
{
  if (!inherits(object, "TriLines"))
    stop("object must be of class \"TriLines\"")
  
  eqn <- object$equation
  xv<-object$x
  yv<-object$y
  nv<-min(6,length(xv),length(yv))
  
  res <- list(txt=object$txt1,
              call=object$call,
              xvec=xv[1:nv],
              yvec=yv[1:nv],
              coefficients=c(object$intercept,object$slope),
              eqn=eqn)
  class(res) <- "summary.TriLines"
  res
}

########################

print.summary.TriLines <- function(x, ...)
{
  
  cat("Call:\n")
  print(x$call)
  cat("\n Selected x points and estimated y points that fall on the Line
      (first 6 or fewer are printed) \n")
  print(x$xvec)
  print(x$yvec)
  
  cat("\nEquation of the", x$txt,"\n")
  Eqn<-x$eqn; names(Eqn)<-c()
  print(Eqn)
  
  cat("\nCoefficients of the line (in the form: y = slope * x + intercept ) \n")
  print(x$coefficients)
}

########################
plot.TriLines<-function (x,xlab="x",ylab="y", ...)
{
  if (!inherits(x, "TriLines"))
    stop("x must be of class \"TriLines\"")
  Tr<-x$tri
  
  A<-Tr[1,]; B<-Tr[2,]; C<-Tr[3,]
  xfence<-abs(A[1]-B[1])*.25 #how far to go at the lower and upper ends in the x-coordinate
  xp<-seq(min(A[1],B[1])-xfence,max(A[1],B[1])+xfence,l=100)
  yp<-x$intercept+x$slope*xp
  
  M<-round(x$cent,2);
  if (x$cent.name=="CC")
  {D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2; #midpoints of the edges
  Ds<-rbind(D1,D2,D3)
  } else
  {Ds<-cp2e.tri(Tr,M)}
  
  Xlim<-range(Tr[,1],xp)
  Ylim<-range(Tr[,2],yp)
  xd<-Xlim[2]-Xlim[1]
  yd<-Ylim[2]-Ylim[1]
  
  plot(Tr,pch=".",main=x$mtxt,xlab=xlab,ylab=ylab,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
  lines(xp,yp,lty=1,col=2)
  polygon(Tr)
  L<-Tr; R<-rbind(M,M,M)
  segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  L<-Ds; R<-rbind(M,M,M)
  segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  
  xp1<-median(xp)#as.numeric(quantile(xp)[2]);
  yp1<-x$intercept+x$slope*xp1
  
  txt<-rbind(Tr,M,Ds,c(xp1,yp1)+c(0,-.1))
  xc<-txt[,1]+c(-.02,.02,.02,.05,.05,-.03,.0,0)
  yc<-txt[,2]+c(.02,.02,.02,.02,0,.02,-.04,0)
  txt.str<-c("A","B","C",x$cent.name,"D1","D2","D3",x$txt2)
  text(xc,yc,txt.str,col=c(rep(1,7),2))
  
}

###############################################
#Auxiliary functions for class Digraphs
###############################################
print.Digraphs <- function(x, ...)
{
  if (!inherits(x, "Digraphs"))
    stop("x must be of class \"Digraphs\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nType:\n")
  print(x$type)
}


########################
summary.Digraphs <- function(object, ...)
{
  if (!inherits(object, "Digraphs"))
    stop("object must be of class \"Digraphs\"")
  
  dimn<-dimension(as.matrix(object$vertices))
  
  typ <- object$type
  xv<-as.matrix(object$vertices)
  yv<-as.matrix(object$tess.points)
  ifelse(!is.null(object$S),sv<-as.matrix(object$S),sv<-NA)
  ifelse(!is.null(object$E), ev<-as.matrix(object$E), ev<-NA)
  nx<-min(6,nrow(xv))
  ny<-min(6,nrow(yv))
  na<-min(6,nrow(sv))
 
  ifelse(!is.na(sv),svec<-sv[1:na,],svec<-NA)
  ifelse(!is.na(ev),evec<-ev[1:na,],evec<-NA)
  
  res <- list(txt=object$txt1,
              call=object$call,
              Vname=object$vert.name, Tname=object$tess.name,
              xvec=xv[1:nx,],
              yvec=yv[1:ny,],
              svec=svec, #ifelse(!is.na(sv),svec=sv[1:na,],svec=NA),
              evec=evec, #ifelse(!is.na(ev),evec=ev[1:na,],evec=NA),
              param=object$parameters,
              type=typ,
              Quant=object$quant,
              dimen=dimn)
  
  class(res) <- "summary.Digraphs"
  res
}
########################

print.summary.Digraphs <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nType of the digraph\n")
  print(x$type)
  
  dimn<-x$dimen
  
  if (dimn==1)
  { cat("\n Vertices of the digraph=",  x$Vname, "\n Partition points of the region=", x$Tname,"\n")
  } else
  {cat("\n Vertices of the digraph (first 6 or fewer are printed) \n")
    print(x$xvec)
    cat("\n Partition points of the region (first 6 or fewer are printed) \n")
    print(x$yvec)}
  
  cat("\n Selected tail (or source) points of the arcs in the digraph
      (first 6 or fewer are printed) \n")
  print(x$svec)
  
  cat("\n Selected head (or end) points of the arcs in the digraph
      (first 6 or fewer are printed) \n")
  print(x$evec)
  
  cat("\nParameters of the digraph\n")
  print(x$param)
  
  cat("\nVarious quantities of the digraph\n")
  print(x$Quant)
}

########################
plot.Digraphs<-function (x, jit=0.1, ...)
{
  if (!inherits(x, "Digraphs"))
    stop("x must be of class \"Digraphs\"")
  
  Xp<-as.matrix(x$vertices)
  Yp<-as.matrix(x$tess.points)
  
  dimn<-dimension(Xp)
  
  nx<-length(Xp); ny<-length(Yp)
  yjit<-runif(nx,-jit,jit)
  
  S<-x$S
  E<-x$E
  if (dimn==1)
  {
    Xlim<-range(Xp,Yp)
    xd<-Xlim[2]-Xlim[1]
    
    plot(cbind(Xp[1],0),main=x$mtxt, xlab=" ", ylab=" ",
         xlim=Xlim+xd*c(-.05,.05),ylim=3*c(-jit,jit),pch=".")
    points(Xp, yjit,pch=".",cex=3)
    abline(h=0,lty=1)
    abline(v=Yp,lty=2)
    if (!is.null(S)) {arrows(S, yjit, E, yjit, length = 0.05, col= 4)}
  } else if (dimn==2 && nrow(Yp)==3)
  {Xlim<-range(Yp[,1],Xp[,1])
  Ylim<-range(Yp[,2],Xp[,2])
  xd<-Xlim[2]-Xlim[1]
  yd<-Ylim[2]-Ylim[1]
  plot(Yp,pch=".",main=x$mtxt,xlab="",ylab="",axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
  polygon(Yp)
  points(Xp,pch=1,col=1)
  if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}
  } else if (dimn==2 && nrow(Yp)>3)
  {
    DT<-tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    
    Xlim<-range(Xp[,1],Yp[,1])
    Ylim<-range(Xp[,2],Yp[,2])
    xd<-Xlim[2]-Xlim[1]
    yd<-Ylim[2]-Ylim[1]
    
    plot(Xp,main=x$mtxt, xlab=" ", ylab=" ",
         xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),pch=".",cex=3)
    plot.tri(DT, add=T, do.points=T)
    if (!is.null(S)) {arrows(S[,1], S[,2], E[,1], E[,2], length = 0.1, col= 4)}
  } else
  {stop('Currently only digraphs with vertices of dimension 1 or 2 are plotted.')}
}

###############################################
#Auxiliary functions for class Patterns
###############################################
print.Patterns <- function(x, ...)
{
  if (!inherits(x, "Patterns"))
    stop("x must be of class \"Patterns\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nType:\n")
  print(x$type)
}


########################
summary.Patterns <- function(object, ...)
{
  if (!inherits(object, "Patterns"))
    stop("object must be of class \"Patterns\"")
  
  typ <- object$type
  xv<-as.matrix(object$gen.points)
  yv<-as.matrix(object$ref.points)
  nx<-min(6,nrow(xv))
  ny<-min(6,nrow(yv))
  Npts<-object$num.points
  
  res <- list(txt=object$desc.pat,
              call=object$call,
              xvec=xv[1:nx,],
              yvec=yv[1:ny,],
              param=object$parameters,
              type=typ,
              num.pts=Npts,
              Xlim=object$xlimit,
              Ylim=object$ylimit,
              pat.desc=object$desc.pat
  )
  
  class(res) <- "summary.Patterns"
  res
}
########################

print.summary.Patterns <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nType of the Pattern\n")
  print(x$type)
  
  cat("\nParameters of the Pattern\n")
  print(x$param)
  
  cat("\nStudy Window\n")
  cat("range in x-coordinate =", x$Xlim,"\n")
  cat("range in y-coordinate =", x$Ylim,"\n")
  
  cat("\n Generated Points from Pattern of",x$pat.desc ," \n (first 6 or fewer are printed) \n")
  print(x$xvec)
  
  cat("\n Reference points (i.e., Y points) in the region=
      (first 6 or fewer are printed) \n")
  print(x$yvec)
  
  cat("\nNumber of points\n nx = number of generated points according to the pattern\n ny = number of reference (i.e. Y) points\n")
  print(x$num.pts)
}

########################
plot.Patterns<-function (x, asp=NA,xlab="x",ylab="y",...)
{
  Y<-x$ref.points
  NY<-nrow(Y)
  Xdt<-x$gen.points
  Xlim<-x$xlimit
  Ylim<-x$ylimit
  
  xf<-.01*(Xlim[2]-Xlim[1])
  yf<-.01*(Ylim[2]-Ylim[1])
  
  plot(Y,asp=asp,pch=16,col=2,lwd=2, xlab=xlab,ylab=ylab,main=x$mtxt,xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf))
  points(Xdt,pch=".",cex=3)
  if (x$tri.Y==TRUE && NY==3)
  {polygon(Y)
  } else if (x$tri.Y==TRUE && NY>3)
  {
    DTY<-tri.mesh(Y[,1],Y[,2],duplicate="remove") #Delaunay triangulation based on Y points
    plot.tri(DTY, add=TRUE, do.points=TRUE,col="blue")
  }
}

###############################################
#Auxiliary functions for class Uniform
###############################################
print.Uniform <- function(x, ...)
{
  if (!inherits(x, "Uniform"))
    stop("x must be of class \"Uniform\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nType:\n")
  print(x$type)
}


########################
summary.Uniform <- function(object, ...)
{
  if (!inherits(object, "Uniform"))
    stop("object must be of class \"Uniform\"")
  
  typ <- object$type
  xv<-as.matrix(object$gen.points)
  yv<-as.matrix(object$tess.points)
  dimn<-dimension(yv)
  nx<-min(6,nrow(xv))
  ny<-min(6,nrow(yv))
  Npts<-object$num.points
  
  res <- list(txt=object$desc.pat,
              call=object$call,
              xvec=xv[1:nx,],
              yvec=yv[1:ny,],
              Support=object$tess.points,
              type=typ,
              num.pts=Npts, txt4points=object$txt4pnts,
              Xlim=object$xlimit,
              Ylim=object$ylimit,
              pat.desc=object$desc.pat,
              dimen=dimn
  )
  
  class(res) <- "summary.Uniform"
  res
}
########################

print.summary.Uniform <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nType of the Pattern\n")
  print(x$type)
  
  cat("\nStudy Window\n")
  cat("range in x-coordinate =", x$Xlim,"\n")
  cat("range in y-coordinate =", x$Ylim,"\n")
  
  Dimn<-x$dimen
  
  if (Dimn==2)
  {
    ny<-x$num.pts[2]
    if (ny<=3)
    {
      cat("\n Vertices of the Support of the Uniform Points\n")
      print(x$yvec)
    } else
    {
      cat("\n Support of the Uniform Points is the convex hull of Ypnts of which first 6 or fewer are printed \n")
      print(x$yvec)
      
      cat("\n The two vectors containing the coordinates of the nodes on the boundary of the convex hull of Ypnts \n")
      Yp<-x$Support
      library(tripack)
      Ydeltri<-tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
      print(convex.hull(Ydeltri) )
    }
  } else if (Dimn==3)
  {
    cat("\n Vertices of the Support of the Uniform Points\n")
    print(x$yvec)
  } else
  {stop('Currently summary.Uniform works for 2D and 3D data')}
  
  cat("\n", x$pat.desc ," \n (first 6 or fewer are printed) \n")
  print(x$xvec)
  
  cat("\nNumber of points\n",x$txt4points,"\n")
  print(x$num.pts)
  
}

########################
plot.Uniform<-function (x,asp=NA,xlab="x",ylab="y",zlab="z", ...)
{
  Y<-x$tess.points
  NY<-nrow(Y)
  Xdat<-x$gen.points
  Xlim<-x$xlimit
  Ylim<-x$ylimit
  Zlim<-x$zlimit
  
  xf<-.01*(Xlim[2]-Xlim[1])
  yf<-.01*(Ylim[2]-Ylim[1])
  zf<-.01*(Zlim[2]-Zlim[1])
  
  yv<-as.matrix(Y)
  Dimn<-dimension(yv)
  
  if (Dimn==2)
  {
    plot(Y,asp=asp,pch=16,col=2,lwd=2, xlab=xlab,ylab=ylab,main=x$mtxt,xlim=Xlim+c(-xf,xf),ylim=Ylim+c(-yf,yf))
    if (!is.null(x$out.region))
    {
      polygon(Y,lty=1,border=2,col=5)
      polygon(x$out.region,lty=2)
      points(Xdat,pch=".",cex=3)
    } else if (NY<=3)
    {
      polygon(Y,lty=1)
      points(Xdat,pch=1)
    } else
    {
      library(tripack)
      DTY<-tri.mesh(Y[,1],Y[,2],duplicate="remove") #Delaunay triangulation based on Y points
      plot.tri(DTY, add=TRUE, do.points=T,pch=16,col="blue")
      points(Xdat,pch=".",cex=3)
    }
  } else if (Dimn==3)
  {
    library(plot3D)
    scatter3D(Xdat[,1],Xdat[,2],Xdat[,3], theta =225, phi = 30,main=x$mtxt,xlab=xlab,ylab=ylab,zlab=zlab, bty = "g",
              xlim=Xlim+c(-xf,xf), ylim=Ylim+c(-yf,yf), zlim=Zlim+c(-zf,zf),pch = 20, cex = 1, ticktype = "detailed")
    #add the vertices of the tetrahedron
    points3D(Y[,1],Y[,2],Y[,3], add = TRUE)
    L<-rbind(Y[1,],Y[1,],Y[1,],Y[2,],Y[2,],Y[3,]); R<-rbind(Y[2,],Y[3,],Y[4,],Y[3,],Y[4,],Y[4,])
    segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
  } else
  {stop('Currently plot.Uniform works for 2D and 3D data')}
  
}

###############################################
#Auxiliary functions for class Extrema
###############################################
print.Extrema <- function(x, ...)
{
  if (!inherits(x, "Extrema"))
    stop("x must be of class \"Extrema\"")
  
  cat("Call:\n")
  print(x$call)
  cat("\nType:\n")
  print(x$type)
}


########################
summary.Extrema <- function(object, ...)
{
  if (!inherits(object, "Extrema"))
    stop("object must be of class \"Extrema\"")
  
  typ <- object$type
  xv<-as.matrix(object$X)
  yv<-as.matrix(object$ROI)
  Extr<-as.matrix(object$Ext)
  nx<-min(6,nrow(xv))
  ny<-min(6,nrow(yv))
  Npts<-object$num.points
  
  res <- list(txt1=object$txt1,
              txt2=object$txt2,
              txt3=object$mtxt,
              call=object$call,
              xvec=xv[1:nx,],
              yvec=yv[1:ny,],
              Support=object$ROI,
              type=typ,
              desc=object$desc,
              extr=Extr,
              d2ref=object$dist2ref,
              Nx=Npts
  )
  
  class(res) <- "summary.Extrema"
  res
}
########################

print.summary.Extrema <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nType of the Extrema\n")
  print(x$type)
  
  cat("\n Extrema Points:",x$desc,"\n")
  print(x$extr)
  cat("\n", x$txt2 ,"\n")
  print(x$d2ref)
  
  cat("\n Vertices of the Support of the Data Points\n")
  print(x$yvec)
  print(x$txt1)
  
  cat("\n Data Points, Dt, of size",x$Nx,":\n (first 6 or fewer are printed) \n")
  print(x$xvec)
  
}

########################
plot.Extrema<-function (x,asp=NA,xlab=" ",ylab=" ",zlab="", ...)
{
  Y<-x$ROI
  Xdt<-x$X
  Cent<-x$cent
  dimn<-dimension(as.matrix(Xdt))
  
  if (dimn==1)
  {
    Xlim<-range(Y,Xdt)
    xd<-Xlim[2]-Xlim[1]
    
    plot(cbind(Xdt,0),main=x$mtxt,xlab=xlab,ylab=ylab,xlim=Xlim+xd*c(-.05,.05))
    abline(h=0)
    abline(v=c(a,b,Mc),col=c(1,1,2),lty=2)
    points(cbind(x$Ext,0),pch=4,col=2)
    text(cbind(c(Y,Cent,x$region.centers),-0.1),c("a","b",x$ncent,x$region.names))
  } else if (dimn==2)
  {
    Xlim<-range(Y[,1],Cent[1],Xdt[,1])
    Ylim<-range(Y[,2],Cent[2],Xdt[,2])
    xd<-Xlim[2]-Xlim[1]
    yd<-Ylim[2]-Ylim[1]
    
    if (!is.null(x$ncent) && x$ncent=="CC") {asp<-1}
    
    plot(Xdt,asp=asp,main=x$mtxt,xlab=xlab,ylab=ylab,axes=TRUE,xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05))
    polygon(Y,lty=1)
    for (i in 1:length(x$regions))
    {polygon(x$regions[[i]],lty=3)}
    points(x$Ext,pty=2,pch=4,col="red")
    
    txt<-rbind(Y,Cent,x$region.centers)
    xc<-txt[,1]
    yc<-txt[,2]
    txt.str<-c("A","B","C",x$ncent,x$region.names)
    text(xc,yc,txt.str)
    
  } else if (dimn==3)
  {
    A<-Y[1,]; B<-Y[2,]; C<-Y[3,]; D<-Y[4,];
    Xlim<-range(Y[,1],Cent[1])
    Ylim<-range(Y[,2],Cent[2])
    Zlim<-range(Y[,3],Cent[3])
    
    xd<-Xlim[2]-Xlim[1]
    yd<-Ylim[2]-Ylim[1]
    zd<-Zlim[2]-Zlim[1]
    
    library(plot3D)
    scatter3D(Xdt[,1],Xdt[,2],Xdt[,3],main=x$mtxt,xlab=xlab,ylab=ylab,zlab="", phi =0,theta=40, bty = "g",
              xlim=Xlim+xd*c(-.05,.05),ylim=Ylim+yd*c(-.05,.05),zlim=Zlim+zd*c(-.05,.05), pch = 20, cex = 1, ticktype = "detailed")
    #add the vertices of the tetrahedron
    points3D(Y[,1],Y[,2],Y[,3], add = TRUE)
    # polygon3D(Y[,1],Y[,2],Y[,3],col=NULL,border=1,add=T)
    L<-rbind(A,A,A,B,B,C); R<-rbind(B,C,D,C,D,D)
    segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lwd=2)
    points3D(x$Ext[,1],x$Ext[,2],x$Ext[,3], pch=4,col="red", add = TRUE)
    
    D1<-(A+B)/2; D2<-(A+C)/2; D3<-(A+D)/2; D4<-(B+C)/2; D5<-(B+D)/2; D6<-(C+D)/2;
    L<-rbind(D1,D2,D3,D4,D5,D6); R<-rbind(Cent,Cent,Cent,Cent,Cent,Cent)
    segments3D(L[,1], L[,2], L[,3], R[,1], R[,2],R[,3], add=TRUE,lty=2)
    
    txt<-rbind(Y,Cent,x$region.centers)
    txt.str<-c("A","B","C","D",x$ncent,x$region.names)
    text3D(txt[,1],txt[,2],txt[,3],txt.str,add=T)
    
  } else
  {stop('plot.Extrema works for 1D, 2D and 3D data')}
  
}
