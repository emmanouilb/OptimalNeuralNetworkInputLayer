#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solutions
  setwd("f:/collocation/chap2");
  source("pde_1a_nlbc.R");
  source("dss044.R");
#
# Grid (in x)
  n=21;xl=0;xu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  D=0.001;ua=298;sig=5.670373e-08;
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=1000;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,2*n);
  for(i in 1:n){
    u0[i  ]=1000;
    u0[i+n]=u0[i];
  }
  ncall=0;
#
# ODE integration
  out=ode(func=pde_1a_nlbc,y=u0,times=tout);
  nrow(out)
  ncol(out)
#
# Arrays for display
  u1=matrix(0,nrow=n,ncol=nout);
  u2=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u1[i,it]=out[it,i+1];
    u2[i,it]=out[it,i+1+n];
  }
  }
#
# Numerical solution
  cat(sprintf("\n      t     x   u1(x,t)   u2(x,t)")); 
  for(it in 1:nout){
  for(i  in 1:n){
    cat(sprintf("\n %6.2f%6.2f%10.2f%10.2f",
      tout[it],x[i],u1[i,it],u2[i,it]));
  }
  }
  matplot(x,u1,type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u1(x,t),u2(x,t)");
  matpoints(x,u2,pch="o",col="black");
#
# Calls to ODE routine
  cat(sprintf("\n\n  ncall = %3d\n",ncall));