  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector  
# of u(x,t)
#  
# ux
  table1=splinefun(x,u);
  ux=table1(x,deriv=1);
#
# uxxx
  uxxx=uxxx7c(xl,xu,n,u);
#
# PDE
  ut=rep(0,n);
  for(i in 1:n){
    ut[i]=-uxxx[i]-6*u[i]*ux[i];
  }
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
}