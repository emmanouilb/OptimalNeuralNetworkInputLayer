#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical, analytical solutions
  setwd("f:/collocation/chap3");
  source("pde_2a.R");
#
# Parameters
  D=1;
#
# Grid (in x,y)
  nx=11;ny=11;xl=0;xu=1;yl=0;yu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(nx-1));
  y=seq(from=yl,to=yu,by=(yu-yl)/(ny-1));
  out=(1-1)*nx+6;
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=0.1;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,nx*ny);
  for(j in 1:ny){
  for(i in 1:nx){
    u0[(j-1)*nx+i]=sin(pi*x[i])*cos(pi*y[j]);
  }
  }
  ncall=0;
#
# ODE integration
  u=u0;t=t0;h=0.0001;
  for(i1 in 1:6){
    ua=sin(pi*x[6])*cos(pi*y[1])*exp(-2*D*(pi^2)*t);
    cat(sprintf("\n t = %6.4f  u(x,y,t) = %8.4f%8.4f",
        t,u[out],ua));
    uplot=matrix(0,nrow=ny,ncol=nx);
    for(j in 1:ny){
    for(i in 1:nx){
      uplot[j,i]=u[(j-1)*nx+i];
    }
    }
    persp(x,y,uplot,theta=45,phi=45,xlim=c(xl,xu),
          ylim=c(yl,yu),zlim=c(-1,1),xlab="x",
          ylab="y",zlab="u(x,y,t)");
  for(i2 in 1:200){
    derv=pde_2a(t,u,parm);
    u=u+derv*h;
    t=t+h;
  }
  }
#
# Calls to ODE routine
  cat(sprintf("\n\n ncall = %3d\n",ncall));