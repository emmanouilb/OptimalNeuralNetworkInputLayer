  pde_1=function(t,u,parms){
#
# Function pde_1 computes the t derivative vector of 
# u(x,t)
#
# BCs
  u[1]=1;
  u[n]=0;
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
  for(i in 1:n){
    f=u[i]*(u[i]-1)*(a-u[i]);
    ut[i]=D*uxx[i]+f;
  }
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }