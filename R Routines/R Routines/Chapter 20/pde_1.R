  pde_1=function(t,u,parm){
#
# Function pde_1 computes the t derivative vectors 
# for u1(x,t), u2(x,t), v1(x,t), v2(x,t)
#  
# One vector to four vectors
  u1=rep(0,n);u2=rep(0,n);
  v1=rep(0,n);v2=rep(0,n);
  for(i in 1:n){
    u1[i]=u[i];
    u2[i]=u[i+n];
    v1[i]=u[i+2*n];
    v2[i]=u[i+3*n];
  } 
#
# u1r
  tabler=splinefun(r,u1);
  u1r=tabler(r,deriv=1);
#
# v1r
  tabler=splinefun(r,v1);
  v1r=tabler(r,deriv=1);
#
# BCs 
  u1r[1]=0;u1r[n]=0;
  v1r[1]=0;v1r[n]=0;
#
# u1rr
  tablerr=splinefun(r,u1r);
  u1rr=tablerr(r,deriv=1);
#
# v1rr
  tablerr=splinefun(r,v1r);
  v1rr=tablerr(r,deriv=1);
#
# PDEs
  u1t=rep(0,n);u2t=rep(0,n);
  v1t=rep(0,n);v2t=rep(0,n);
  for(i in 1:n){
    u1t[i]=u2[i];
    v1t[i]=v2[i];
    term1=-exp(-2*u1[i])*(v2[i]^2-v1r[i]^2);
    term2=-2*(u1r[i]*v1r[i]-u2[i]*v2[i]);
    if(i==1){
       u2t[i]=2*u1rr[i]+c1*term1;
       v2t[i]=2*v1rr[i]+c2*term2;}
    if(i> 1){
       u2t[i]=u1rr[i]+(1/r[i])*u1r[i]+c1*term1;
       v2t[i]=v1rr[i]+(1/r[i])*v1r[i]+c2*term2;}
  }
#
# Four vectors to one vector
  ut=rep(0,4*n);
  for(i in 1:n){
    ut[i]  =u1t[i];
    ut[i+n]=u2t[i];
    ut[i+2*n]=v1t[i];
    ut[i+3*n]=v2t[i];
  }  
#
# Increment calls to pde_1
  ncall <<- ncall+1;
#
# Return derivative vector
  return(list(c(ut)));
  }