#
# Boussinesq
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap11");
  source("pde_1a.R");source("ua_1.R");
  source("uat_1.R");
#
# Grid in x
  n=101;xl=-15;xu=35;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  c=0.9;
#
# Independent variable for ODE integration
  nout=4;t0=0;tf=15;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# ICs from analytical solution
  u0=rep(0,2*n);
  for(i in 1:n){
    u0[i]  = ua_1(x[i],t0);
    u0[i+n]=uat_1(x[i],t0);
  }
  ncall=0;
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
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
# Arrays for analytical solution, errors in numerical solution 
  u1a=matrix(0,nrow=n,ncol=nout);
  err=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
    u1[1,it]=0;
    u1[n,it]=0;
    for(i in 1:n){
       u1a[i,it]=ua_1(x[i],t[it]);       
       err[i,it]=u1[i,it]-u1a[i,it];
    }     
  }
#
# Display selected output
  cat(sprintf(" xl = %4.1f  xu= %4.1f  c = %4.1f\n",xl,xu,c));
  for(it in 1:nout){
    cat(sprintf("\n     t      x   u1(x,t)  u1a(x,t)  err(x,t)\n"));
    iv=seq(from=1,to=n,by=10);
    for(i in iv){
      cat(sprintf("%6.2f%7.1f%10.6f%10.6f%10.6f\n",
              t[it],x[i],u1[i,it],u1a[i,it],err[i,it]));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Plot 2D numerical solution
    matplot(x,u1,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u(x,t)",main="Boussinesq");
    matpoints(x,u1a,pch="o",col="black");
#
# Plot 3D numerical solution
    persp(x,t,u1,theta=55,phi=45,xlim=c(xl,xu),
          ylim=c(t0,t[nout]),xlab="x",ylab="t",
          zlab="u(x,t)");