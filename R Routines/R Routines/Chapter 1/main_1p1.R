#
# Previous workspaces are cleared
  rm(list=ls(all=TRUE))
#
# Define uniform grid
  xl=0;xu=1;n=11;
  x=seq(from=xl,to=xu,by=(xu-xl)/(n-1));
#
# Define function to be approximated
  u=sin(pi*x);
#
# Set up spline table
  utable=splinefun(x,u);
#
# Compute spline approximation
  us=utable(x);
#
# Display comparison of function and its spline
# approximation
  cat(sprintf("\n    x         u        us        diff"));
  for(i in 1:n){
    cat(sprintf("\n%5.2f%10.5f%10.5f%12.7f",
      x[i],u[i],us[i],us[i]-u[i]));
  }