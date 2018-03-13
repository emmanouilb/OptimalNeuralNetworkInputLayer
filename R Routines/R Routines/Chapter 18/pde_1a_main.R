#
# Euler-Poisson-Darboux
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator, routines
  library(deSolve);
  setwd("f:/collocation/chap18");
  source("pde_1a.R");
#
# Select case
  ncase=1;
#
# Parameters
  if(ncase==1){lam=0;c=5;};
  if(ncase==2){lam=1;c=5;};
#
# Spatial interval
  n=51;r0=1;
  r=seq(from=0,to=r0,by=(r0-0)/(n-1));
#
# Time interval
  t0=0;tf=0.5;nout=6;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
  ncall=0;
#
# Initial conditions (ICs)
  u0=rep(0,2*n);
  for(i in 1:n){
    u0[i]=exp(-c*r[i]^2);
    u0[i+n]=0;
  }
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Store numerical solution for plotting
  u1=matrix(0,nrow=n,ncol=nout);
  u2=matrix(0,nrow=n,ncol=nout);
   t=rep(0,nout);
  for(it in 1:nout){
    t[it]=out[it,1];
  for(i in 1:n){
    u1[i,it]=out[it,i+1];
    u2[i,it]=out[it,i+1+n];
  }
  }
#
# Numerical solution
  for(it in 1:nout){
  cat(sprintf("\n      t       r   u(r,t)")); 
  iv=seq(from=1,to=n,by=5);
  for(i in iv){
    cat(sprintf("\n %6.3f%8.3f%9.4f",t[it],r[i],u1[i,it]));
  }
    cat(sprintf("\n"));
  }
  cat(sprintf("\n ncall = %4d\n",ncall));
#
# Plot numerical solution
  matplot(r,u1,type="l",lwd=2,col="black",lty=1,
    xlab="r",ylab="u1(r,t)",main="Euler-Poisson-Darboux");