  pde_1=function(t,u,parm){
#
# Function pde_1 computes the t derivative vectors 
# for u1(x,t), u2(x,t)
#  
# One vector to two vectors
  u1=rep(0,n);u2=rep(0,n);
  for(i  in 1:n){
    u1[i]=u[i];
    u2[i]=u[i+n];
  } 
#
# ux
  tablex=splinefun(x,u1);
  u1x=tablex(x,deriv=1);
#
# BCs 
  u1x[1]=0;u1x[n]=0;
#
# uxx
  tablexx=splinefun(x,u1x);
  u1xx=tablexx(x,deriv=1);
#
# PDEs
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=u2[i];
    u2t[i]=-a*u1xx[i]-b*u1[i];
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