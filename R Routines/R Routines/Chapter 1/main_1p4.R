#
# Previous workspaces are cleared
  rm(list=ls(all=TRUE))
#
# Define grid
  xl=0;xu=1;n=11;ng=2;
#
# Uniform grid
  if(ng==1){
    x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
    for(i in 2:n){
      cat(sprintf("\n i = %2d  x = %6.4f  dx = %6.4f",
                 i,x[i],x[i]-x[i-1]));
    }
  }
#
# Variable grid
  if(ng==2){
    x=rep(0,n);
    x[1]=0;xin=(1/0.2020)/(n-1);
    cat(sprintf("\n x[1] = %6.4f  xin = %6.4f",x[1],xin));
    for(i in 2:n){
      x[i]=x[i-1]+xin/i;
      cat(sprintf("\n i = %2d  x = %6.4f  dx = %6.4f",
                  i,x[i],x[i]-x[i-1]));
    } 
  }
#
# Define function to be approximated, and its
# derivatives
     a=2;
     u=exp(a*x);
    ux=a*exp(a*x);
   uxx=a^2*exp(a*x);
  uxxx=a^3*exp(a*x);
#
# Set up spline table for function, derivatives
     utable=splinefun(x,u);
#
# Compute spline approximation of function
# derivatives
      us=utable(x);
     usx=utable(x,deriv=1);
    usxx=utable(x,deriv=2);
   usxxx=utable(x,deriv=3);
#
# Display comparison of function and its spline
# approximation, and its derivatives
#
# u
  cat(sprintf("\n      x         u        us        diff"));
  for(i in 1:n){
    cat(sprintf("\n%7.4f%10.5f%10.5f%12.7f",
      x[i],u[i],us[i],us[i]-u[i]));
  }
#
# ux
  cat(sprintf("\n      x        ux       usx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%7.4f%10.5f%10.5f%12.7f",
      x[i],ux[i],usx[i],usx[i]-ux[i]));
  }
#
# uxx
  cat(sprintf("\n      x       uxx      usxx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%7.4f%10.5f%10.5f%12.7f",
      x[i],uxx[i],usxx[i],usxx[i]-uxx[i]));
  }
#
# uxxx
  cat(sprintf("\n      x      uxxx     usxxx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%7.4f%10.5f%10.5f%12.7f",
      x[i],uxxx[i],usxxx[i],usxxx[i]-uxxx[i]));
  }