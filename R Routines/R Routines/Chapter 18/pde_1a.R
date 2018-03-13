  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vectors
# for u1(r,t), u2(r,t)
#
# One vector to two vectors
  u1=rep(0,n);u2=rep(0,n);
  for(i in 1:n){
   u1[i]=u[i];
   u2[i]=u[i+n];
   }
#
# u1r
  table1=splinefun(r,u1);
  u1r=table1(r,deriv=1);
#
# Neumann BCs
  u1r[1]=0;u1r[n]=0;
#
# u1rr
  table2=splinefun(r,u1r);
  u1rr=table2(r,deriv=1);
#
# PDE
  u1t=rep(0,n);u2t=rep(0,n);
  for(i in 1:n){
    if(t< 1.0e-04){
      if(i==1){
        u1t[i]=u2[i];
        u2t[i]=1/(1+lam)*3*u1rr[i];
      }
      if(i> 1){
        u1t[i]=u2[i];
        u2t[i]=1/(1+lam)*(u1rr[i]+(2/r[i])*u1r[i]);
      }
    }
    if(t>=1.0e-04){
      if(i==1){
        u1t[i]=u2[i];
        u2t[i]=-lam/t*u2[i]+3*u1rr[i];
      }
      if(i> 1){
        u1t[i]=u2[i];
        u2t[i]=-lam/t*u2[i]+(u1rr[i]+(2/r[i])*u1r[i]);
      }
    }
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
}