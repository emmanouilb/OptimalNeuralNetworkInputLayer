  ua_1=function(x,t){
#
# Function ua_1 computes the analytical solution of 
# the Boussinesq equation for comparison with the 
# numerical solution
#
# Analytical solution
  xi=(1-c^2)^(0.5)/2*(x-c*t);
  ua=(3/2)*(1-c^2)*(cosh(xi))^(-2);
  return(c(ua));
  }