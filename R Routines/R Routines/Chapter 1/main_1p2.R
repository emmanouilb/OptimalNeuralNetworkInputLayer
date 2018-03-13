#
# Previous workspaces are cleared
  rm(list=ls(all=TRUE))
#
# Define uniform grid
  xl=0;xu=1;n=11;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Define function to be approximated, and its
# derivatives
     u=sin(pi*x);
    ux=pi*cos(pi*x);
   uxx=-pi^2*sin(pi*x);
  uxxx=-pi^3*cos(pi*x);
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
  cat(sprintf("\n    x         u        us        diff"));
  for(i in 1:n){
    cat(sprintf("\n%5.2f%10.5f%10.5f%12.7f",
      x[i],u[i],us[i],us[i]-u[i]));
  }
#
# ux
  cat(sprintf("\n    x        ux       usx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%5.2f%10.5f%10.5f%12.7f",
      x[i],ux[i],usx[i],usx[i]-ux[i]));
  }
#
# uxx
  cat(sprintf("\n    x       uxx      usxx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%5.2f%10.5f%10.5f%12.7f",
      x[i],uxx[i],usxx[i],usxx[i]-uxx[i]));
  }
#
# uxxx
  cat(sprintf("\n    x      uxxx     usxxx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%5.2f%10.5f%10.5f%12.7f",
      x[i],uxxx[i],usxxx[i],usxxx[i]-uxxx[i]));
  }