  pde_1a=function(t,u,parm){
#
# Function pde_1a computes the t derivative vector for 
# a PDE with a mixed partial derivative
#
# BCs
  u[1]=sin(pi*x[1]/(xu-xl))*exp(a*t);
  u[n]=sin(pi*x[n]/(xu-xl))*exp(a*t);
#
# PDE
  ut=rep(0,n);
  ut=solve(cm)%*%u;
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
  }