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
     a=1;
     u=exp(a*x);
   uxx=a^2*exp(a*x);
#
# Set up spline tables for u, ux; compute uxx
# by stagewise differentiation
  utable=splinefun(x,u);
  usx=utable(x,deriv=1);
  uxtable=splinefun(x,usx);
  usxx=uxtable(x,deriv=1);
#
# Display comparison of exact and spline derivatives
#
# uxx
  cat(sprintf("\n     x       uxx      usxx        diff"));
  for(i in 1:n){
    cat(sprintf("\n%6.4f%10.5f%10.5f%12.7f",
      x[i],uxx[i],usxx[i],usxx[i]-uxx[i]));
  }