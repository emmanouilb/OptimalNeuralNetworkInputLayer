  pde_1=function(t,u,parms){
#
# Function pde_1 computes the t derivative vector 
# of u(x,t)
#
# BCs
  u[1]=ua_1(x[1],t);
  u[n]=ua_1(x[n],t);
#
# ux
  tablex=splinefun(x,u);
  ux=tablex(x,deriv=1);
#
# uxx
  tablexx=splinefun(x,ux);
  uxx=tablexx(x,deriv=1);
#
# PDE
  ut=rep(0,n);
  for(i in 2:(n-1)){
    ut[i]=-(u[i]^2)*ux[i]+uxx[i]+(2/3)*u[i]^3*(1-u[i]^2);
  }
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }