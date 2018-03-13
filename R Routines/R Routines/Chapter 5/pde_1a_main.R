#
# Korteweg-deVries
#
# Delete previous workspaces
  rm(list=ls(all=TRUE))
#
# Access ODE integrator
  library("deSolve");
#
# Access functions for numerical solution
  setwd("f:/collocation/chap5");
  source("pde_1a.R");source("inital_1.R");
  source("uxxx7c.R");source("simp.R")    ;
  source("ua.R")    ;source("dss004.R")  ;
#
# Select case; 1 - one soliton; 2 - two solitons;
  ncase=1;
  if(ncase==1){c1=1;n=201;}
  if(ncase==2){c1=2;c2=0.5;n=301;}
#
# Spatial grid
  xl=-30;xu=70;dx=(xu-xl)/(n-1);
  x=seq(from=xl,to=xu,by=dx);
#
# Independent variable for ODE integration
  t0=0;tf=30;nout=4;
  tout=seq(from=t0,to=tf,by=(tf-t0)/(nout-1));
  ncall=0;
#
# Initial condition 
  u0=inital_1(x,t0);
#
# ODE integration
  if(ncase==1){
    out=lsodes(func=pde_1a,y=u0,times=tout,
               sparsetype ="sparseint");}
  if(ncase==2){
    out=lsodes(func=pde_1a,y=u0,times=tout,
               sparsetype ="sparseint",
               rtol=1e-10,atol=1e-10,maxord=5);}
  nrow(out)
  ncol(out)
#
# Output, ncase=1   
  if(ncase==1){  
# 
#   Store analytical solution, errors in numerical solution
    u=matrix(0,nrow=n,ncol=nout);  
  u_a=matrix(0,nrow=n,ncol=nout); 
  err=matrix(0,nrow=n,ncol=nout); 
    for(it in 1:nout){
      for(i in 1:n){
          u[i,it]=out[it,i+1];
        u_a[i,it]=ua(x[i],tout[it]);       
        err[i,it]=u[i,it]-u_a[i,it];
      }    
    }
#
#   Display selected output
    cat(sprintf("\n ncase = %2d  c1 = %5.2f\n",ncase,c1));
    for(it in 1:nout){
      cat(sprintf("\n     t       x        u(it,i)      u_a(it,i)      
                  err(it,i)\n"));
      iv=seq(from=1,to=n,by=5);
      for(i in iv){  
        cat(sprintf("%6.2f%8.3f%15.6f%15.6f%15.6f\n",
              tout[it],x[i],u[i,it],u_a[i,it],err[i,it]));
      }
#
#     Calculate and display three invariants
      ui=u[,it];  
      uint=simp(xl,xu,n,ui);  
      cat(sprintf("\n Invariants at t = %5.2f",tout[it]));
      cat(sprintf("\n    I1 = %10.4f   Mass conservation"    ,
                  uint[1]));
      cat(sprintf("\n    I2 = %10.4f   Energy conservation"  ,
                  uint[2]));
      cat(sprintf("\n    I3 = %10.4f   Whitham invariant\n\n",
                  uint[3]));
    }
    cat(sprintf("    ncall = %4d\n\n",ncall));
#
#   Plot numerical and analytical solutions
    matplot(x,u,type="l",lwd=2,col="black",lty=1,
            xlab="x",ylab="u(x,t)",
            main="Korteweg-deVries");
    matpoints(x,u_a,pch="o",col="black");
#
# End of ncase=1
  }  
#
# Output, ncase=2  
  if(ncase==2){
#       
#   Store numerical solution
    u=matrix(0,nrow=n,ncol=nout); 
    for(it in 1:nout){
      for(i in 1:n){
        u[i,it]=out[it,i+1];
      }    
    }
#
#   Display selected output
    cat(sprintf("\n ncase = %2d   c1 = %5.2f   c2 = %5.2f\n",
                ncase,c1,c2));
    for(it in 1:nout){
      cat(sprintf("     t         x     u(i,it)\n"));
      iv=seq(from=1,to=n,by=5);
      for(i in iv){ 
        cat(sprintf("%6.2f%10.3f%12.6f\n",tout[it],x[i],u[i,it]));
      }
      cat(sprintf("\n"));
#
#     Calculate and display three invariants
      ui=u[,it];  
      uint=simp(xl,xu,n,ui);  
      cat(sprintf("\n Invariants at t = %5.2f",tout[it]));
      cat(sprintf("\n    I1 = %10.4f   Mass conservation"    ,
                  uint[1]));
      cat(sprintf("\n    I2 = %10.4f   Energy conservation"  ,
                  uint[2]));
      cat(sprintf("\n    I3 = %10.4f   Whitham invariant\n\n",
                  uint[3]));
    }
    cat(sprintf("  ncall = %4d\n\n",ncall));
#
#   Plot numerical solution
    par(mfrow=c(2,2))
    for(it in 1:nout){
      if(it==1){
        plot(x,u[,1],type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u(x,t=0)",main="u(x,t=0");}
      if(it==2){
        plot(x,u[,2],type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u(x,t=10)",main="u(x,t=10)");}
      if(it==3){
        plot(x,u[,3],type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u(x,t=20)",main="u(x,t=20)");}
      if(it==4){
        plot(x,u[,4],type="l",lwd=2,col="black",lty=1,
          xlab="x",ylab="u(x,t=30)",main="u(x,t=30)");}
    }
#
# End of ncase=2
  }