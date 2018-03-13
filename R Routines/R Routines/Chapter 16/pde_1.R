  pde_1=function(t,u,parm){
#
# Function pde_1 computes the t derivative vector 
# for u1(x,t) u2(x,t) 
#  
# One vector to two vectors
  u1=rep(0,n);u2=rep(0,n);
  for(i in 1:n){
    u1[i]=u[i];
    u2[i]=u[i+n];
  } 
#
# u1x
  tablex=splinefun(x,u1);
  u1x=tablex(x,deriv=1);
#
# u2x
  tablex=splinefun(x,u2);
  u2x=tablex(x,deriv=1);
#
# BCs 
  u1x[1]=0;u1x[n]=0;
  u2x[1]=0;u2x[n]=0;
#
# Nonlinear term
  u1u2x=rep(0,n);
  for(i in 1:n){
    u1u2x[i]=2*u2[i]/u1[i]*u1x[i];
  }
#
# u1u2xx
  tablex=splinefun(x,u1u2x);
  u1u2xx=tablex(x,deriv=1);
#
# u2xx
  tablexx=splinefun(x,u2x);
  u2xx=tablexx(x,deriv=1);
#
# PDEs
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=-k*u2[i];
    u2t[i]=D*(u2xx[i]-u1u2xx[i]);  
  }
#
# Two vectors to one vector
  ut=rep(0,2*n);
  for(i in 1:n){
    ut[i]  =u1t[i];
    ut[i+n]=u2t[i];
  }  
#
# Increment calls to pde_1
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
  }