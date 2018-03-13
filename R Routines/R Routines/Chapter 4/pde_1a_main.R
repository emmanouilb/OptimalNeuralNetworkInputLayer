#
# Burgers
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap4");
  source("pde_1a.R");source("phi.R");
#
# Parameters
  xl=0;xu=1;
  vis=0.003;
#
# Spatial grid
  xl=0;xu=1;n=201;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1)); 
#
# Independent variable for ODE integration
  t0=0;tf=1;nout=11;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
  ncall=0;
#
# Initial condition 
  u0=rep(0,n);
  for(i in 1:n){
    u0[i]=phi(x[i],t0);
  }
#
# ODE integration
  out=lsodes(func=pde_1a,y=u0,times=tout,
             sparsetype ="sparseint")
  nrow(out)
  ncol(out)
#
# Arrays for display
   u=matrix(0,nrow=n,ncol=nout);
  ua=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u[i,it]=out[it,i+1];
   ua[i,it]=phi(x[i],tout[it]);
  }
  }
#
# Numerical, analytical solutions
  cat(sprintf("\n      t     x    u(x,t)   ua(x,t)        diff")); 
  for(it in 1:nout){
  iv=seq(from=1,to=n,by=10);
  for(i  in iv){
    cat(sprintf("\n %6.2f%6.2f%10.4f%10.4f%12.6f",
      tout[it],x[i],u[i,it],ua[i,it],u[i,it]-ua[i,it]));
  }
  }
  matplot(x,u,type="l",lwd=2,col="black",lty=1,
    xlab="x",ylab="u(x,t)",main="Burgers");
  matpoints(x,ua,pch="o",col="black");
#
# Calls to ODE routine
  cat(sprintf("\n\n  ncall = %3d\n",ncall));