#
# Function pde_1c computes the t derivative vector of 
# u1(x,t),u2(x,t)
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
# u1xx
  tablexx=splinefun(x,u1x);
  u1xx=tablexx(x,deriv=1);
#
# u2xx
  tablexx=splinefun(x,u2x);
  u2xx=tablexx(x,deriv=1);
#
# PDEs
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=-v1*u1x[i]+D1*u1xx[i]-k1*u1[i]+k2*u2[i];
    u2t[i]=-v2*u2x[i]+D2*u1xx[i]+k1*u1[i]-k2*u2[i];
  }
#
# Increment calls to pde_1c
  ncall <<- ncall+1; 
#
# Two vectors to one vector
  ut=rep(0,2*n);
  for(i in 1:n){
    ut[i]  =u1t[i];
    ut[i+n]=u2t[i];
  }
#
# Return derivative vector
  return(list(c(ut)));
  }