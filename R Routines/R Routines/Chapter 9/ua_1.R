  ua_1=function(x,t){
#
# Function ua_1 computes the exact solution of the 
# Fisher-Kolmogorov equation for comparison with the 
# numerical solution 
#
# Analytical solution
  z=x-c*t;
  ua=1/(1+a*exp(b*z))^s;
#
# Return solution
  return(c(ua));
  }