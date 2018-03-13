  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector 
# of u(x,t)
#
# ux
  tablex=splinefun(x,u);
  ux=tablex(x,deriv=1);
#
# BCs
  ux[1]=0;ux[n]=0;
#
# uxx
  tablexx=splinefun(x,ux);
  uxx=tablexx(x,deriv=1);
#
# PDE
  ut=rep(0,n);
  for(i in 1:n){
    ut[i]=-(u[i]^2)*ux[i]+uxx[i]+(2/3)*u[i]^3*(1-u[i]^2);
  }
#
# Increment calls to pde_1a
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }