  pde_1a=function(t,u,parm){
#
# Function pde_1a computes the t derivative vectors
# for u1(x,t), u2(x,t)
#  
# One vector to two vectors
  u1=rep(0,n);u2=rep(0,n);
  for(i  in 1:n){
    u1[i]=u[i];
    u2[i]=u[i+n];
  } 
#
# BCs 
  u1[1]=ua_1(x[1],t);
  u1[n]=ua_1(x[n],t);
#
# uxx
  tablexx=splinefun(x,u1);
  u1xx=tablexx(x,deriv=2);
#
# PDE
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=u2[i];
    u2t[i]=-a*u1xx[i]-b*u1[i]-g*u1[i]^m;
  }
  u1t[1]=0;u1t[n]=0;
#
# Two vectors to one vector
  ut=rep(0,2*n);
  for(i in 1:n){
    ut[i]  =u1t[i];
    ut[i+n]=u2t[i];
  }  
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
  }