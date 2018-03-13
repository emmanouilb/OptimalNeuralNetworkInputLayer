  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector of 
# u(x,t)
#
# BCs
  u[1]=ua_1(x[1],t);
  u[n]=ua_1(x[n],t);
#
# uxx
  tablexx=splinefun(x,u);
  uxx=tablexx(x,deriv=2);
#
# PDE
  ut=rep(0,n);
  for(i in 1:n){
    ut[i]=uxx[i]+u[i]*(1-u[i]^q);
  }
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1a
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }