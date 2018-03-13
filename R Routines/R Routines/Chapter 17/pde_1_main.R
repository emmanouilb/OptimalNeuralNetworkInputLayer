#
# Fitzhugh-Nagumo
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap17");
  source("pde_1.R");
#
# Grid in x
  n=101;xl=-60;xu=20;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  ncase=1;
  if(ncase==1){a=1;D=1  ;}
  if(ncase==2){a=1;D=0.1;}
#
# Independent variable for ODE integration
  nout=4;t0=0;tf=60;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,n);sr=1/sqrt(2*D);
  for(i in 1:n){
    u0[i]=1/(1+exp(sr*x[i]));
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
# Arrays for numerical solution
   u=matrix(0,nrow=n,ncol=nout);
  ua=matrix(0,nrow=n,ncol=nout);
  t=rep(0,nout);
  for(it in 1:nout){
  for(i  in 1:n){
      t[it]=out[it,1];   
    u[i,it]=out[it,i+1];
   ua[i,it]=1/(1+exp(sr*x[i]+(a-0.5)*t[it]));
  }     
  }
#
# Display selected output
  cat(sprintf("\n a = %4.2f,  D = %4.2f\n",a,D));
  for(it in 1:nout){
    cat(sprintf("\n     t        x      u(x,t)     ua(x,t)   diff\n"));
    u[1,it]=1;u[n,it]=0;
    iv=seq(from=1,to=n,by=10);
    for(i in iv){
      diff=u[i,it]-ua[i,it];
      cat(sprintf("%6.2f%9.3f%12.6f%12.6f%12.6f\n",
              t[it],x[i],u[i,it],ua[i,it],diff));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Graphical output
#
# u, ua
    matplot(x,u,type="l",lwd=2,col="black",lty=1,
      xlab="x",ylab="u(x,t)",main="Fitzhugh-Nagumo");
    matpoints(x,ua,pch="o",col="black");