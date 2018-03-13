#
# Fisher-Kolmogorov
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap9");
  source("pde_1a.R");source("ua_1.R");
#
# Grid in x
  n=31;xl=-5;xu=10;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Parameters
  q=1;
  a=2^(1/2)-1;
  b=q/(2*(q+2))^(1/2);
  c=(q+4)/(2*(q+2))^(1/2);
  s=2/q; 
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=5;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# Initial condition
  u0=rep(0,n);
  for(i in 1:n){
    u0[i]=ua_1(x[i],t0);
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
# Arrays for analytical solution, errors 
# in numerical solution 
   ua=matrix(0,nrow=n,ncol=nout);
  err=matrix(0,nrow=n,ncol=nout);
  for(it in 1:nout){
    u[1,it]=ua_1(x[1],t[it]);
    u[n,it]=ua_1(x[n],t[it]);
    for(i in 1:n){
       ua[i,it]=ua_1(x[i],t[it]);       
      err[i,it]=u[i,it]-ua[i,it];
    }     
  }
#
# Display selected output
  cat(sprintf("\n q = %4.2f,  a = %4.2f,  b = %4.2f,  
              c = %4.2f  s = %4.2f\n",q,a,b,c,s));
  for(it in 1:nout){
    cat(sprintf("\n     t       x      u(x,t)     ua(x,t)    
                err(x,t)\n"));
    iv=seq(from=1,to=n,by=5);
    for(i in iv){
      cat(sprintf("%6.2f%8.3f%12.6f%12.6f%12.6f\n",
                  t[it],x[i],u[i,it],ua[i,it],err[i,it]));
    }
    cat(sprintf("\n")); 
  }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Graphical output
  matplot(x,u,type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u(x,t)",main="Fisher-Kolmogorov");
  matpoints(x,ua,pch="o",col="black");