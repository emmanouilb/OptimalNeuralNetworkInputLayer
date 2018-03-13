#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap3");
  source("pde_3a_dnr.R");
#
# Parameters
  D=1;ua=1;
#
# Grid (in x,y)
  nx=11;ny=11;nz=11;
  xl=0;xu=1;yl=0;yu=1;zl=0;zu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(nx-1));
  y=seq(from=yl,to=yu,by=(yu-yl)/(ny-1));
  z=seq(from=zl,to=zu,by=(zu-zl)/(nz-1));
  out=(6-1)*nx*ny+(6-1)*nx+6;
#
# Independent variable for ODE integration
  nout=21;t0=0;tf=0.5;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,nx*ny*nz);
  ncall=0;
#
# ODE integration
   tplot=rep(0,nout);
   uplot=rep(0,nout);
  u=u0;t=t0;h=0.0005;
  for(i1 in 1:nout){
    cat(sprintf("\n t = %6.4f  u(0.5,0.5,0.5,t) = %8.4f",
        t,u[out]));
    tplot[i1]=t;
    uplot[i1]=u[out];
  for(i2 in 1:50){
    derv=pde_3a_dnr(t,u,parm);
    u=u+derv*h;
    t=t+h;
  }
  }
#
# Plot u(x=0.5,y=0.5,z=0.5,t) against t
  par(mfrow=c(1,1))
  plot(tplot,uplot,
    xlab="t",ylab="u(0.5,0.5,0.5,t)",
    type="l",lwd=2);
#
# Calls to ODE routine
  cat(sprintf("\n\n ncall = %3d\n",ncall));