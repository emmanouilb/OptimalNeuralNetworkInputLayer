  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vectors
# of u1(x,t),u2(x,t)
#
# One vector to two vectors
  u1=rep(0,n);u2=rep(0,n);
  for(i in 1:n){
   u1[i]=u[i];
   u2[i]=u[i+n];
   }
#
# u1x
  table1=splinefun(x,u1);
  u1x=table1(x,deriv=1);
#
# Neumann BCs
  u1x[1]=0;u1x[n]=0;
#
# u1xx
  table2=splinefun(x,u1x);
  u1xx=table2(x,deriv=1);
#
# PDE
  mueps=1/(mu*eps);
  sigeps=sigma/eps;
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=u2[i];
    u2t[i]=mueps*u1xx[i]-sigeps*u2[i];
  }
#
# Two vectors to one vector
  ut=rep(0,2*n);
  for(i in 1:n){  
      ut[i]=u1t[i];
    ut[i+n]=u2t[i];
  }  
#
# Increment calls to pde_1a
  ncall <<- ncall+1;
#
# Return derivative vector            
  return(list(c(ut)))
#
# End of pde_1a
}