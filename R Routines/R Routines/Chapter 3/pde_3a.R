  pde_3a=function(t,u1d,parms){
#
# Function pde_3a computes the t derivative vector 
# of u(x,y,z,t)
#
# 1D vector to 3D array
  u=array(u1d,c(nx,ny,nz));
#
# ux
  ux=array(0,c(nx,ny,nz));
  for(k in 1:nz){
  for(j in 1:ny){
    tablex=splinefun(x,u[k,j,]);
    ux[k,j,]=tablex(x,deriv=1);
  }
  }
#
# uy
  uy=array(0,c(nx,ny,nz));
  for(k in 1:nz){
  for(i in 1:nx){
    tabley=splinefun(y,u[k,,i]);
    uy[k,,i]=tabley(y,deriv=1);
  }
  }
#
# uz
  uz=array(0,c(nx,ny,nz));
  for(j in 1:ny){
  for(i in 1:nx){
    tablez=splinefun(z,u[,j,i]);
    uz[,j,i]=tablez(y,deriv=1);
  }
  }
#
# BCs, x=0,xL
  ux[,, 1]=0;
  ux[,,nx]=0;
#
# BCs, y=0,yL
  uy[, 1,]=0;
  uy[,ny,]=0;
#
# BCs, z=0,zL
  uz[ 1,,]=0;
  uz[nz,,]=0;
#
# uxx
  uxx=array(0,c(nx,ny,nz));
  for(k in 1:nz){
  for(j in 1:ny){
    tablexx=splinefun(x,ux[k,j,]);
    uxx[k,j,]=tablexx(x,deriv=1);
  }
  }
#
# uyy
  uyy=array(0,c(nx,ny,nz));
  for(k in 1:nz){
  for(i in 1:nx){
    tableyy=splinefun(y,uy[k,,i]);
    uyy[k,,i]=tableyy(y,deriv=1);
  }
  }
#
# uzz
  uzz=array(0,c(nx,ny,nz));
  for(j in 1:ny){
  for(i in 1:nx){
    tablezz=splinefun(z,uz[,j,i]);
    uzz[,j,i]=tablezz(y,deriv=1);
  }
  }
#
# PDE
  ut=array(0,c(nx,ny,nz));
  ut=uxx+uyy+uzz;
#
# 3D array to 1D vector
   u1dt=ut[,,];
#
# Increment calls to pde_3a
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(c(u1dt));
  }