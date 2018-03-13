#
# Cahn-Hilliard
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap12");
  source("pde_1a.R");
#
# Select case (for IC)
  ncase=3;
#
# Grid in x
  n=51;xl=-5;xu=5;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  D=1;gam=0.5;
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=5;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# IC
  u0=rep(0,n);
  for(i in 1:n){
    if(ncase==1){u0[i]=0};
    if(ncase==2){u0[i]=tanh(x[i]/sqrt(2*gam))};
    if(ncase==3){
      ul=-1;uu=1;
      if(i<11){u0[i]=ul};
      if(i>41){u0[i]=uu};
      if((i>=11)&(i<=41)){
        u0[i]=ul+(uu-ul)*(i-11)/(41-11);
#       cat(sprintf("\n %5d%10.4f%10.4f",i,x[i],u0[i]));
      }
    }
  }
  ncall=0;
#
# ODE integration
  out=lsodes(y=u0,times=tout,func=pde_1a,
      sparsetype ="sparseint",rtol=1e-6,atol=1e-6,
      maxord=5);
  nrow(out)
  ncol(out)
#
# Arrays for numerical solution
  u=matrix(0,nrow=n,ncol=nout);
  t=rep(0,nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u[i,it]=out[it,i+1];
      t[it]=out[it,1];       
  } 
  }
#
# Display selected output
  cat(sprintf(" xl = %4.1f  xu= %4.1f\n  D = %4.1f  gam = %4.1f\n",
              xl,xu,D,gam));
  for(it in 1:nout){
    cat(sprintf("\n      t      x    u(x,t)\n"));
    iv=seq(from=1,to=n,by=5);
    for(i in iv){
      cat(sprintf("%7.3f%7.3f%10.6f\n",t[it],x[i],u[i,it]));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Equilibrium solution
  if(ncase>1){
    ue=rep(0,n);
    cat(sprintf("\n      Equiilbrium solution\n"));
    for(i in 1:n){
      ue[i]=tanh(x[i]/sqrt(2*gam));
    }
    iv=seq(from=1,to=n,by=5);
    for(i in iv){
      cat(sprintf("\n       %7.3f%10.6f",x[i],ue[i]));
    }
  }
#
# Plot 2D numerical solution
    matplot(x,u,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u(x,t)",main="Cahn-Hilliard");
    if(ncase>1){
    matpoints(x,ue,pch="o",col="black");
    }
#
# Plot 3D numerical solution
    persp(x,t,u,theta=0,phi=55,xlim=c(xl,xu),
          ylim=c(t0,t[nout]),xlab="x",ylab="t",