  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the derivative vector in t
#  
# ux
  table1=splinefun(x,u);
  ux=table1(x,deriv=1);
#
# BC at x = 0
  ux[1]=0;
#
# BC at x = 1
  ux[n]=0;
#
# uxx
  table2=splinefun(x,ux);
  uxx=table2(x,deriv=1);
#
# PDE
  ut=rep(0,n);
  for(i in 1:n){
    ut[i]=vis*uxx[i]-u[i]*ux[i];
  }
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
}