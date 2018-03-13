#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap3");
  source("pde_3a.R");
#
# Parameters
  D=1;
#
# Grid (in x,y)
  nx=11;ny=11;nz=11;
  xl=0;xu=1;yl=0;yu=1;zl=0;zu=1;
  x=seq(from=xl,to=xu,by=(xu-xl)/(nx-1));
  y=seq(from=yl,to=yu,by=(yu-yl)/(ny-1));
  z=seq(from=zl,to=zu,by=(zu-zl)/(nz-1));
  out=(1-1)*nx*ny+(1-1)*nx+1;
#
# Independent variable for ODE integration
  nout=21;t0=0;tf=0.1;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,nx*ny*nz);
  for(k in 1:nz){
  for(j in 1:ny){
  for(i in 1:nx){
    u0[(k-1)*nx*ny+(j-1)*nx+i]=cos(pi*x[i])*
                               cos(pi*y[j])*
                               cos(pi*z[k]);
  }
  }
  }
  ncall=0;
#
# ODE integration
   tplot=rep(0,nout);
   uplot=rep(0,nout);
  uaplot=rep(0,nout);
  u=u0;t=t0;h=0.0001;
  for(i1 in 1:nout){
    ua=cos(pi*x[1])*cos(pi*y[1])*cos(pi*z[1])*exp(-3*D*(pi^2)*t);
    cat(sprintf("\n t = %6.4f  u(0,0,0,t) = %8.4f%8.4f",
        t,u[out],ua));
    tplot[i1]=t;
    uplot[i1]=u[out];
   uaplot[i1]=ua;
  for(i2 in 1:50){
    derv=pde_3a(t,u,parm);
    u=u+derv*h;
    t=t+h;
  }
  }
#
# Plot u(x=0,y=0,z=0,t) against t
  par(mfrow=c(1,1))
  plot(tplot,uplot,xlab="t",ylab="u(0,0,0,t)",
       type="l",lwd=2);
  points(tplot,uaplot,pch="o",lwd=2);
#
# Calls to ODE routine
  cat(sprintf("\n\n ncall = %3d\n",ncall));