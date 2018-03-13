  pde_1a_nbc=function(t,u,parms){
#
# Function pde_1a_nbc computes the t derivative vectors 
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
  ut=D*uxx;
#
# Increment calls to pde_1a_nbc
  ncall <<- ncall+1; 
#
# Return derivative vector
  return(list(c(ut)));
  }