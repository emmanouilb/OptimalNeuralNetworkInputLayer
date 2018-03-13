#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap2");
  source("pde_1a.R");
#
# Step through ncase
  for(ncase in 1:2){
#
# Grid (in x)
  n=21;xl=0;xu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  D=1;
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=0.25;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  if(ncase==1)u0=sin(pi*x);
  if(ncase==2)u0=rep(1,n);
  ncall=0;
#
# ODE integration
  out=ode(func=pde_1a,y=u0,times=tout);
#
# Arrays for display
   u=matrix(0,nrow=n,ncol=nout);
  ua=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u[i,it]=out[it,i+1];
    if(ncase==1){ua[i,it]=sin(pi*x[i])*exp(-D*(pi^2)*
                 tout[it]);}
  }
  }
  cat(sprintf("\n ncase = %2d",ncase));
#
# Numerical, analytical solutions, ncase = 1
  if(ncase==1){
    cat(sprintf("\n      t     x    u(x,t)   ua(x,t)        
                 diff")); 
    for(it in 1:nout){
    for(i  in 1:n){
      cat(sprintf("\n %6.2f%6.2f%10.4f%10.4f%12.6f",
        tout[it],x[i],u[i,it],ua[i,it],u[i,it]-ua[i,it]));
    }
    }
    matplot(x,u,type="l",lwd=2,col="black",lty=1,
            xlab="x",ylab="u(x,t)");
    matpoints(x,ua,pch="o",col="black");
#
# End ncase=1
  }
#
# Numerical solution, ncase = 2
  if(ncase==2){
  cat(sprintf("\n      t     x    u(x,t)")); 
    for(it in 1:nout){
    for(i  in 1:n){
      cat(sprintf("\n %6.2f%6.2f%10.4f",
        tout[it],x[i],u[i,it]));
 
    }
    }
    matplot(x[2:(n-1)],u[2:(n-1),],type="l",lwd=2,col="black",
            lty=1,xlab="x",ylab="u(x,t)");
#
# End ncase=2
  }
#
# Calls to ODE routine
  cat(sprintf("\n\n  ncall = %3d\n",ncall));
#
# Next case
  }