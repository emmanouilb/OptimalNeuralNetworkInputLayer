#
# Maxwell
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator, routines
  library(deSolve);
  setwd("f:/collocation/chap6");
  source("pde_1a.R");
#
# Select case
#
# ncase=1 - no damping; ncase=2 - damping
  ncase=1;
#
# Parameters
  if(ncase==1){eps=1;mu=1;sigma=0;}
  if(ncase==2){eps=1;mu=1;sigma=1;}
  Re_lam=-0.5*sigma/eps;
  Im_lam= 0.5*sqrt(4*pi^2/(mu*eps)-(sigma/eps)^2);
#
# Spatial interval
  n=51;xl=0;xu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Time interval
  t0=0;tf=1;nout=11;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
  ncall=0;
#
# Initial conditions (ICs)
  u0=rep(0,2*n);
  u0[1:n]=cos(pi*x);
  u0[(n+1):(2*n)]=Re_lam*u0[1:n];
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Store numerical, analytical solutions for plotting
   u1=matrix(0,nrow=n,ncol=nout);
   u2=matrix(0,nrow=n,ncol=nout);
  u1a=matrix(0,nrow=n,ncol=nout);
  t=rep(0,nout);
  for(it in 1:nout){
    t[it]=out[it,1];
  for(i  in 1:n){
    u1[i,it]=out[it,i+1];
    u2[i,it]=out[it,i+1+n];
   u1a[i,it]=exp(Re_lam*tout[it])*cos(Im_lam*tout[it])*
             cos(pi*x[i]);
  }
  }
#
# Numerical, analytical solutions
  cat(sprintf("\n      t     x    u(x,t)   ua(x,t)        
              diff")); 
  for(it in 1:nout){
  iv=seq(from=1,to=n,by=5);
  for(i  in iv){
    cat(sprintf("\n %6.2f%6.2f%10.4f%10.4f%12.6f",
      tout[it],x[i],u1[i,it],u1a[i,it],u1[i,it]-u1a[i,it]));
  }
  }
  cat(sprintf("\n ncall = %4d\n",ncall));
#
# Plot numerical, analytical solutions
  matplot(x,u1,type="l",lwd=2,col="black",lty=1,
    xlab="x",ylab="u(x,t)",main="Maxwell");
  matpoints(x,u1a,pch="o",col="black");
#
# Plot 3D numerical solution
  persp(x,tout,u1,theta=55,phi=45,xlim=c(xl,xu),
        ylim=c(t0,tout[nout]),zlim=c(-1.1,1.1),xlab="x",
        ylab="t",zlab="u(x,t)");