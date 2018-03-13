  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector
# for u(x,t)
#
# Dirichlet BCs at x = xl,xu
  u[1]=ua_1(x[1],t);
  u[n]=ua_1(x[n],t);
#
# ux
  tablex=splinefun(x,u);
  ux=tablex(x,deriv=1);
#
# Neumann BCs at x = xl,xu
  ux[1]=uax_1(x[1],t);
  ux[n]=uax_1(x[n],t);
#
# uxx
  tablexx=splinefun(x,ux);
  uxx=tablexx(x,deriv=1);
#
# uxxx
  tablexxx=splinefun(x,uxx);
  uxxx=tablexxx(x,deriv=1);
#
# uxxxx
  tablexxxx=splinefun(x,uxxx);
  uxxxx=tablexxxx(x,deriv=1);
#
# PDE
  ut=rep(0,n);
  for(i in 2:(n-1)){
    ut[i]=-u[i]*ux[i]-alpha*uxx[i]-
           beta*uxxx[i]-gamma*uxxxx[i];
  }
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector            
  return(list(c(ut)))
}