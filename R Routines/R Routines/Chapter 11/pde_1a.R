  pde_1a=function(t,u,parm){
#
# Function pde_1a computes the t derivative vector 
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
  u1[1]=0;
  u1[n]=0;
#
# uxx
  tablexx=splinefun(x,u1);
  u1xx=tablexx(x,deriv=2);
#
# BCs 
  u1xx[1]=0;
  u1xx[n]=0;
#
# uxxxx
  tablexxxx=splinefun(x,u1xx);
  u1xxxx=tablexxxx(x,deriv=2);
#
# u1^2
  u1s=rep(0,n);
  for(i in 1:n){
    u1s[i]=u1[i]^2;
  }
#
# u1sxx
  tablexx=splinefun(x,u1s);
  u1sxx=tablexx(x,deriv=2);
#
# PDE
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=u2[i];
    u2t[i]=u1xx[i]-u1xxxx[i]-u1sxx[i];
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