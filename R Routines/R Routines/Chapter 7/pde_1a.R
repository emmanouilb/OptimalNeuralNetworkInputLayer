  pde_1a=function(t,u,parms){
#
# Function pde_1a computes the t derivative vector of 
# c(x,t), phi(x,t), ci(t)
#
# One vector to two vectors and a scalar
  c=rep(0,n);phi=rep(0,n);
  for(i in 1:n){
    c[i]=u[i];
  phi[i]=u[i+n];
  }
  ci=u[2*n+1];
#
# BCs
  c[1]=1;phi[1]=1;
#
# cx
  table=splinefun(x,c);
  cx=table(x,deriv=1);
#
# BC
  cx[n]=km*(ci-c[n]);
#
# cxx
  table=splinefun(x,cx);
  cxx=table(x,deriv=1);
#
# phix
  table=splinefun(x,phi);
  phix=table(x,deriv=1);
#
# BC
  phix[n]=0;
#
# phixx
  table=splinefun(x,phix);
  phixx=table(x,deriv=1);
#
# c*phix
  cphix=rep(0,n);
  for(i in 1:n){
    cphix[i]=c[i]*phix[i];
  }
#
# (c*phix)x
  table=splinefun(x,cphix);
  cphixx=table(x,deriv=1);
#
# PDEs
  ct=rep(0,n);phit=rep(0,n);
  for(i in 1:n){
      ct[i]=cxx[i]-c11*cx[i]+c12*cphixx[i];
    phit[i]=(1/mu)*(phixx[i]-c21*c[i]);
  }
  ct[1]=0;phit[1]=0;
  cit=(1/tau)*(c[n]-ci);
#
# Increment calls to pde_1a
  ncall <<- ncall+1; 
#
# Two vectors and a scalar to one vector
  ut=rep(0,2*n+1);
  for(i in 1:n){
      ut[i]=ct[i];
    ut[i+n]=phit[i];
  }
  ut[2*n+1]=cit;
#
# Return derivative vector
  return(list(c(ut)));
  }