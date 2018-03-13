  pde_1=function(t,u,parm){
#
# Function pde_1 computes the t derivative vector 
# for u(x,t)
#
# ux
  tablex=splinefun(x,u);
  ux=tablex(x,deriv=1);
#
# uxx
  tablexx=splinefun(x,u);
  uxx=tablexx(x,deriv=2);
#
# uxxx
  tablexxx=splinefun(x,ux);
  uxxx=tablexxx(x,deriv=2);
#
# RHS group
  ug=rep(0,n);
  for(i in 1:n){
    ug[i]=-2*k*ux[i]-3*u[i]*ux[i]+
           2*ux[i]*uxx[i]+u[i]*uxxx[i];
  }
#
# PDE
  ut=rep(0,n);
  ut=solve(cm)%*%ug;
#
# Increment calls to pde_1
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
  }