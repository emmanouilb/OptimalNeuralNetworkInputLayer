#
# Poisson-Nernst-Planck
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap7");
  source("pde_1a.R");
#
# Parameters
     u=0.001;   xL=1.0e-08;     D=1.0e-10;
         Z=1; e=1.6022e-19; kB=1.3806e-23;
       T=300;   c0=1.0e+02;    E0=5.0e+06;
      mu=0.1;   phi0=5e-02; eps=710.0e-12;
        km=1;     tau=0.01;   
#
# Constants
  c11=u*xL/D;
  c12=((Z*e)/(kB*T))*(xL*c0*E0);
  c21=((Z*e)/eps)*((xL^2*c0)/phi0);
#
# Select case
  ncase=2;
#
# No electrodiffusion
  if(ncase==1){
    n=21;c12=0;
    t0=0;tf=0.1
  }
#
# Electrodiffusion
  if(ncase==2){
    n=21;
    t0=0;tf=0.2
  }
  cat(sprintf("\n c11 = %10.3e\n c12 = %10.3e\n c21 = %10.3e\n",
              c11,c12,c21));
#
# Grid (in x)
  xl=0;xu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Independent variable for ODE integration
  nout=6;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial conditions
  c=rep(0,n);phi=rep(0,n);
  c[1]=1;phi[1]=1;
  u0=rep(0,2*n+1);
  for(i in 1:n){
   u0[i]=c[i];
   u0[i+n]=phi[i];
  }
  u0[2*n+1]=0;
  ncall=0;
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Arrays for display
    c=matrix(0,nrow=n,ncol=nout);
  phi=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
  for(i  in 1:n){
      c[i,it]=out[it,i+1];
    phi[i,it]=out[it,i+1+n];
  }
  }
#
# Numerical solutions
  cat(sprintf("\n          t         x    c(x,t)  phi(x,t)")); 
    for(it in 1:nout){
    for(i  in 1:n){
      cat(sprintf("\n %10.5f%10.4f%10.4f%10.4f",
        tout[it],x[i],c[i,it],phi[i,it]));
 
    }
      cat(sprintf("\n"));
    }
    matplot(x,c,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="c(x,t)",main="Poisson-Nernst-Planck");
    matplot(x,phi,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="phi(x,t)",main="Poisson-Nernst-Planck");
#
# Calls to ODE routine
  cat(sprintf("\n\n  ncall = %3d\n",ncall));