#
# Keller-Segel
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap16");
  source("pde_1.R");source("u1a.R");
  source("u2a.R");
#
# Grid in x
  n=101;xl=-10;xu=15;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  k=1;D=1;c=1;
  cat(sprintf("\n\n k = %5.2f   D = %5.2f   c = %5.2f\n",
              k,D,c));
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=5;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition from analytical solutions (t=0)
  u0=rep(0,2*n);
  for(i in 1:n){
    u0[i]  =u1a(x[i],0);
    u0[i+n]=u2a(x[i],0);
  } 
  ncall=0;
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Arrays for plotting numerical, analytical solutions
  u1=matrix(0,nrow=n,ncol=nout);
  u2=matrix(0,nrow=n,ncol=nout);
 ua1=matrix(0,nrow=n,ncol=nout);
 ua2=matrix(0,nrow=n,ncol=nout);
   t=rep(0,nout);
  for(it in 1:nout){
    t[it]=out[it,1];
    for(i in 1:n){
       u1[i,it]=out[it,i+1];
       u2[i,it]=out[it,i+1+n];
      ua1[i,it]=u1a(x[i],t[it]);
      ua2[i,it]=u2a(x[i],t[it]);
    }
  }
#
# Display selected output
#
# Step through t
  for(it in 1:nout){
    cat(sprintf("\n     t       x   u1(x,t)    ua1(x,t)   diff1\n"));
    cat(sprintf(  "     t       x   u2(x,t)    ua2(x,t)   diff2\n"));
#
#   Step through x
    iv=seq(from=1,to=n,by=10);
    for(i in iv){
      diff1=u1[i,it]-ua1[i,it];
      cat(sprintf("%6.2f%8.3f%10.6f%12.6f%12.6f\n",
              t[it],x[i],u1[i,it],ua1[i,it],diff1));
      diff2=u2[i,it]-ua2[i,it];
      cat(sprintf("%6.2f%8.3f%10.6f%12.6f%12.6f\n",
              t[it],x[i],u2[i,it],ua2[i,it],diff2));
    cat(sprintf("\n")); 
    }
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Plot 2D numerical solutions
    matplot(x,u1,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u1(x,t)",main="Keller-Segel");
    matpoints(x,ua1,pch="o",col="black");
    matplot(x,u2,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u2(x,t)",main="Keller-Segel");
    matpoints(x,ua2,pch="o",col="black");