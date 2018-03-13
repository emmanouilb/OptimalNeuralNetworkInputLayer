#
# Camassa-Holm
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap13");
  source("pde_1.R");
#
# Grid in x
  n=119;xl=-5;xu=7;dx=(xu-xl)/(n+1);
  x=seq(from=(xl+dx),to=(xu-dx),by=dx);
#
# Parameters
  c=0.5;k=1;
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=0.5;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# IC
  u0=rep(0,n);
  for(i in 1:n){
    u0[i]=exp(-c*x[i]^2);
  }
  ncall=0;
#
# Coefficient matrix
  cm=matrix(0,nrow=n,ncol=n);
  for(i in 1:n){
  for(j in 1:n){
    if(i==j){cm[i,j]=(1+2/dx^2)};
    if(abs(i-j)==1){cm[i,j]=-1/dx^2};
    if(abs(i-j)>1){cm[i,j]=0}; 
  }
  }
#
# ODE integration
# out=lsodes(y=u0,times=tout,func=pde_1,
#     sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
#     maxord=5);
  out=ode(y=u0,times=tout,func=pde_1);
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
# Display selected output
  cat(sprintf("   xl = %4.1f   xu = %4.1f\
    c = %4.1f    k = %4.1f\n",xl,xu,c,k));
  for(it in 1:nout){
    cat(sprintf("\n      t       x    u(x,t)\n"));
    iv=seq(from=1,to=n,by=10);
    for(i in iv){
      cat(sprintf("%7.2f%8.2f%10.5f\n",
                  t[it],x[i],u[i,it]));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Arrays for plotting
  up=matrix(0,nrow=(n+2),ncol=nout);
  xp=rep(0,(n+2));xp[1]=-5;xp[n+2]=7;
  for(i  in 1:n){xp[i+1]=x[i]};
  for(it in 1:nout){
    up[1,it]=0;up[n+2,it]=0;
    for(i in 1:n){up[i+1,it]=u[i,it]};
  }
#
# Plot numerical solution
  matplot(xp,up,type="l",lwd=2,col="black",lty=1,
    xlab="x",ylab="u(x,t)",main="Camassa-Holm");