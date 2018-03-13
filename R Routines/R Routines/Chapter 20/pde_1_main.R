#
# Einstein-Maxwell
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap20");
  source("pde_1.R");
#
# Grid in r
  n=51;rl=0;ru=1;
  r=seq(from=rl,to=ru,by=(ru-rl)/(n-1));
#
# Parameters
  ncase=1;
  if(ncase==1){c=10;c1=0;c2=0;}
  if(ncase==2){c=10;c1=1;c2=1;}
#
# Independent variable for ODE integration
  nout=6;t0=0;tf=1;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
#
# ICs 
  u0=rep(0,4*n);
  for(i in 1:n){
    u0[i]    =exp(-c*r[i]^2);
    u0[i+n]  =0;
    u0[i+2*n]=exp(-c*r[i]^2);
    u0[i+3*n]=0;
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
  u1=matrix(0,nrow=n,ncol=nout);
  u2=matrix(0,nrow=n,ncol=nout);
  v1=matrix(0,nrow=n,ncol=nout);
  v2=matrix(0,nrow=n,ncol=nout);
  t=rep(0,nout);
  for(it in 1:nout){
  for(i  in 1:n){
    u1[i,it]=out[it,i+1];
    u2[i,it]=out[it,i+1+n];
    v1[i,it]=out[it,i+1+2*n];
    v2[i,it]=out[it,i+1+3*n];
       t[it]=out[it,1];       
  }     
  }
#
# Display selected output
  cat(sprintf("\n   c = %4.2f,  c1 = %4.2f,  c2 = %4.2f\n",
              c,c1,c2));
  for(it in 1:nout){
    iv=seq(from=1,to=n,by=10);
    for(i in iv){
    cat(sprintf("\n     t       r     u1(r,t)     u2(r,t)"));
    cat(sprintf("\n     t       r     v1(r,t)     v2(r,t)\n"));
      cat(sprintf("%6.2f%8.3f%12.6f%12.6f\n",
              t[it],r[i],u1[i,it],u2[i,it]));
      cat(sprintf("%6.2f%8.3f%12.6f%12.6f\n",
              t[it],r[i],v1[i,it],v2[i,it]));
    }
   }  
  cat(sprintf(" ncall = %4d\n",ncall));
#
# Plot 2D numerical solution
    matplot(r,u1,type="l",lwd=2,col="black",lty=1,
      xlab="r",ylab="u1(r,t)",main="Einstein-Maxwell");
    matplot(r,u2,type="l",lwd=2,col="black",lty=1,
      xlab="r",ylab="u2(r,t)",main="Einstein-Maxwell");
    matplot(r,v1,type="l",lwd=2,col="black",lty=1,
      xlab="r",ylab="v1(r,t)",main="Einstein-Maxwell");
    matplot(r,v2,type="l",lwd=2,col="black",lty=1,
      xlab="r",ylab="v2(r,t)",main="Einstein-Maxwell");