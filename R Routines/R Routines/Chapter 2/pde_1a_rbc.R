  pde_1a_rbc=function(t,u,parms){
#
# Function pde_1a_rbc computes the t derivative vectors 
# of u(x,t)
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
# BCs
  u2x=rep(0,n);
  u1x[1]=u1[1]-ua;u1x[n]=ua-u1[n]; 
  u2x[1]=u2[1]-ua;u2x[n]=ua-u2[n]; 
#
# u1xx
  tablexx=splinefun(x,u1x);
  u1xx=tablexx(x,deriv=1);
#
# u2xx
  nl=2;nu=2;
  u2xx=dss044(xl,xu,n,u2,u2x,nl,nu);
#
# PDE
  u1t=D*u1xx;
  u2t=D*u2xx;
#
# Two vectors to one vector
  ut=rep(0,2*n);
  for(i in 1:n){
      ut[i]=u1t[i];
    ut[i+n]=u2t[i];
  }
#
# Increment calls to pde_1a_rbc
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }