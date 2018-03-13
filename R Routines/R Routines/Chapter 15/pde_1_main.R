#
# Gierer-Meinhardt
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap15");
  source("pde_1.R");
#
# Grid in x
  n=101;xl=0;xu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  ncase=1;
  D1=1;D2=1;c=10;
  mu=1;nu=1;
  if(ncase==1){r=1 ;}
  if(ncase==2){r=10;}
#
# Independent variable for ODE integration
  nout=11;t0=0;tf=0.1;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# ICs
  u0=rep(0,2*n);
  for(i in 1:n){
    u0[i]  =exp(-c*x[i]^2);
    u0[i+n]=exp(-c*x[i]^2);
  }
  ncall=0;
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Arrays for numerical solution
  u1=matrix(0,nrow=n,ncol=nout);
  u2=matrix(0,nrow=n,ncol=nout);
  t=rep(0,nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u1[i,it]=out[it,i+1];
    u2[i,it]=out[it,i+1+n];
       t[it]=out[it,1];       
  }     
  }
#
# Display selected output
  cat(sprintf("\n r = %3.1f  mu = %3.1f  nu = %3.1f\n",r,mu,nu));
  for(it in 1:nout){
    cat(sprintf("\n     t        x   u1(x,t)   u2(x,t)\n"));
    iv=seq(from=1,to=n,by=10);
    for(i in iv){
      cat(sprintf("%6.2f%9.3f%10.6f%10.6f\n",
              t[it],x[i],u1[i,it],u2[i,it]));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Plot 2D numerical solution
    matplot(x,u1,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u1(x,t)",main="Gierer-Meinhardt");
    matplot(x,u2,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u2(x,t)",main="Gierer-Meinhardt");