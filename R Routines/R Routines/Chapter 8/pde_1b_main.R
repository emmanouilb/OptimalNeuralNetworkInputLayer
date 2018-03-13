#
# Fokker-Planck
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap8");
  source("pde_1b.R");
#
# Grid (in x)
  n=151;xl=-5;xu=10;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  D=0.1;v=1;c=2;
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=5;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,n);
  for(i in 1:n){
    u0[i]=exp(-c*x[i]^2);
  }
  ncall=0;
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1b,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Arrays for display
   u=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u[i,it]=out[it,i+1];
  }
  }
#
# Numerical solution
  cat(sprintf("\n      t     x    u(x,t)")); 
  for(it in 1:nout){
  iv=seq(from=1,to=n,by=10);
  for(i  in iv){
    cat(sprintf("\n %6.2f%6.2f%10.4f",
        tout[it],x[i],u[i,it]));
    }
    cat(sprintf("\n")); 
    }
#
# Graphical output
  matplot(x,u,type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u(x,t)",main="Fokker-Planck");
#
# Calls to ODE routine
  cat(sprintf("\n\n  ncall = %3d\n",ncall));