  pde_2a=function(t,u1d,parms){
#
# Function pde_2a computes the t derivative vector
# of u(x,y,t)
#
# 1D vector to 2D array
  u=matrix(0,nrow=ny,ncol=nx);
  for(j in 1:ny){
  for(i in 1:nx){
    u[j,i]=u1d[(j-1)*nx+i];  
  }
  }
#
# BCs, x=0,xL
  u[, 1]=0;
  u[,nx]=0;
#
# uy
  uy=matrix(0,nrow=ny,ncol=nx);
  for(i in 1:nx){
    tabley=splinefun(y,u[,i]);
    uy[,i]=tabley(y,deriv=1);
  }  
#
# BCs, y=0,yL
  uy[ 1,]=0;
  uy[ny,]=0;
#
# uxx
  uxx=matrix(0,nrow=ny,ncol=nx);
  for(j in 1:ny){
    tablexx=splinefun(x,u[j,]);
    uxx[j,]=tablexx(x,deriv=2);
  }
#
# uyy
  uyy=matrix(0,nrow=ny,ncol=nx);
  for(i in 1:nx){
    tableyy=splinefun(y,uy[,i]);
    uyy[,i]=tableyy(y,deriv=1);
  }
#
# PDE
  ut=matrix(0,nrow=ny,ncol=nx);
  ut=uxx+uyy;
#
# BCs, x=0,xL
  ut[, 1]=0;
  ut[,nx]=0;
#
# 2D array to 1D vector
  u1dt=rep(0,nx*ny);
  for(j in 1:ny){
  for(i in 1:nx){
    u1dt[(j-1)*nx+i]=ut[j,i];  
  }
  }
#
# Increment calls to pde_2a
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(c(u1dt));
  }