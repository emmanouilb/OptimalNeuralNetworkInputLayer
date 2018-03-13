  pde_1a=function(t,u,parm){
#
# Function pde_1a computes the t derivative vector 
# for u(x,t)
#
# BCs 
  u[1]=-1;
  u[n]= 1;
#
# uxx
  tablexx=splinefun(x,u);
  uxx=tablexx(x,deriv=2);
#
# BCs 
  uxx[1]=0;
  uxx[n]=0;
#
# uxxxx
  tablexxxx=splinefun(x,uxx);
  uxxxx=tablexxxx(x,deriv=2);
#
# u^3-u-gam*uxx
  ug=rep(0,n);
  for(i in 1:n){
    ug[i]=u[i]^3-u[i]-gam*uxx[i]
  }
#
# (u^3-u-gam*uxx)xx
  tablegxx=splinefun(x,ug);
  ugxx=tablegxx(x,deriv=2);
#
# PDE
  ut=rep(0,n);
  for(i in 2:(n-1)){
    ut[i]=D*ugxx[i];
  }
  ut[1]=0;ut[n]=0;
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
  }