#
# Kuramoto-Sivashinsky
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator, routines
  library(deSolve);
  setwd("C:/Users/eurico/SunspotAnalysis/SpatialTemporalFeatureSelectionOptimalApproach/R Routines/R Routines/Chapter 19");
  source("pde_1a.R");source("ua_1.R");
  source("uax_1.R");
#
# Select case
  ncase=1;
  if(ncase==1){
    alpha=1;gamma=1;c0=1;
    #beta=4*(alpha*gamma)^0.5;
    beta=1;
    k=(alpha/gamma)^0.5;
    lambda=-c0*k-(3/2)*beta*k^3;
    c1=(15/76)*(16*alpha-beta^2/gamma)+15*beta*k+60*gamma*k^2;
    c2=-(15*beta+180*gamma*k);
    c3=60*gamma;
#
#   Spatial grid
    n=200;xl=0;xu=1;
    x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
#   Independent variable for ODE integration
    t0=0;tf=0.15;nout=6;
    tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
    ncall=0;    
  } 
#
  if(ncase==2){
    alpha=1;beta=0;gamma=1;c0=0;
    k=(11/19)^0.5;
#
#   Spatial grid
    n=101;xl=-10;xu=20;
    x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
#   Independent variable for ODE integration
    t0=0;tf=10;nout=6;
    tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
    ncall=0;    
  }
#
# Initial condition (IC)
  u0=rep(0,n);
  for(i in 1:n){
    u0[i]=ua_1(x[i],0);
  }
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Store numerical, analytical solutions for plotting
   u=matrix(0,nrow=n,ncol=nout);
  ua=matrix(0,nrow=n,ncol=nout);
   t=rep(0,nout);
  for(it in 1:nout){
    t[it]=out[it,1];
  for(i in 1:n){
    u[i,it]=out[it,i+1];
   ua[i,it]=ua_1(x[i],t[it]);
  }
   u[1,it]=ua_1(x[1],t[it]);
   u[n,it]=ua_1(x[n],t[it]);
  }
#
# Numerical solution
  for(it in 1:nout){
  cat(sprintf("\n      t       x   u(x,t)  ua(x,t)     diff")); 
  if(ncase==1){iv=seq(from=1,to=n,by=5);}
  if(ncase==2){iv=seq(from=1,to=n,by=10);}
  for(i in iv){
    diff=u[i,it]-ua[i,it];
    cat(sprintf("\n %6.3f%8.3f%9.4f%9.4f%9.4f",
                t[it],x[i],u[i,it],ua[i,it],diff));
  }
    cat(sprintf("\n"));
  }
  cat(sprintf("\n ncall = %4d\n",ncall));
#
# Plot numerical solution
  matplot(x,u,type="l",lwd=2,col="black",lty=1,
    xlab="x",ylab="u(x,t)",main="Kuramoto-Sivashinsky");
    matpoints(x,ua,pch="o",col="black");
#
# Plot 3D numerical solution
    persp(x,t,u,theta=55,phi=45,xlim=c(xl,xu),
          ylim=c(t0,t[nout]),xlab="x",ylab="t",
          zlab="u(x,t)");