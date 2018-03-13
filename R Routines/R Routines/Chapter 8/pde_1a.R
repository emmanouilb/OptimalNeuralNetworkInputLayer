  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector of
# u(x,t)
#
# BCs
  u[1]=0;u[n]=0;
#
# ux
  tablex=splinefun(x,u);
  ux=tablex(x,deriv=1);
#
# uxx
  tablexx=splinefun(x,u);
  uxx=tablexx(x,deriv=2);
#
# PDE
  ut=rep(0,n);
  for(i in 2:(n-1)){
    ut[i]=-v*ux[i]+D*uxx[i];
  }
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1a
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }