  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector of 
# u(x,t)
#
# BCs
  u[1]=0;u[n]=0; 
#
# uxx
  table=splinefun(x,u);
  uxx=table(x,deriv=2);
#
# PDE
  ut=D*uxx;
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }