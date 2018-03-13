#
# Test PDE with mixed partial derivative
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap13");
  source("pde_1a.R");
#
# Grid in x
  n=39;xl=0;xu=1;dx=(xu-xl)/(n+1);
  x=seq(from=(xl+dx),to=(xu-dx),by=dx);
#
# Parameters
  a=1/(1-(pi/(xu-xl))^2);
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=20;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# IC
  u0=rep(0,n);
  for(i in 1:n){
    u0[i]=sin(pi*x[i]/(xu-xl));
  }
  ncall=0;
#
# Coefficient matrix
  cm=matrix(0,nrow=n,ncol=n);
  for(i in 1:n){
  for(j in 1:n){
    if(i==j){cm[i,j]=(1-2/dx^2)};
    if(abs(i-j)==1){cm[i,j]=1/dx^2};
    if(abs(i-j)>1){cm[i,j]=0}; 
  }
  }
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Arrays for numerical solution
  u=matrix(0,nrow=n,ncol=nout);
  t=rep(0,nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u[i,it]=out[it,i+1];
      t[it]=out[it,1];       
  } 
  }
#       
# Arrays for analytical solution, errors in numerical solution 
   ua=matrix(0,nrow=n,ncol=nout);
  err=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
    u[1,it]=sin(pi*x[1]/(xu-xl))*exp(a*t[it]);
    u[n,it]=sin(pi*x[n]/(xu-xl))*exp(a*t[it]);
    for(i in 1:n){
      ua[i,it]=sin(pi*x[i]/(xu-xl))*exp(a*t[it]);    
     err[i,it]=u[i,it]-ua[i,it];
    }     
  }
#
# Display selected output
  cat(sprintf(" xl = %4.1f   xu = %4.1f  a = %4.1f\n",xl,xu,a));
  for(it in 1:nout){
    cat(sprintf("\n     t      x    u(x,t)   ua(x,t)  err(x,t)\n"));
    iv=seq(from=1,to=n,by=1);
    for(i in iv){
      cat(sprintf("%6.1f%7.3f%10.6f%10.6f%10.6f\n",
              t[it],x[i],u[i,it],ua[i,it],err[i,it]));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Arrays for plotting
   up=matrix(0,nrow=(n+2),ncol=nout);
  uap=matrix(0,nrow=(n+2),ncol=nout);
   xp=rep(0,(n+2));xp[1]=0;xp[n+2]=1;
   for(i in 1:n){xp[i+1]=x[i]};
  for(it in 1:nout){
    up[1,it]=0;up[n+2,it]=0;
    for(i in 1:n){up[i+1,it]=u[i,it]};
    uap[1,it]=0;uap[n+2,it]=0;
    for(i in 1:n){uap[i+1,it]=ua[i,it]};
  }
#
# Plot 2D numerical, analytical solutions
    matplot(xp,up,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u(x,t)",main="Mixed partial");
    matpoints(xp,uap,pch="o",col="black");