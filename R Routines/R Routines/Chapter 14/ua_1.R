  ua_1=function(x,t){
#
# Function ua_1 computes the exact solution of the 
# Burgers-Huxley equation for comparison with the 
# numerical solution 
#
# Analytical solution
  expp=exp( (1/3)*x+(1/9)*t);
  expm=exp(-(1/3)*x-(1/9)*t);
  ua=((1/2)*(1+(expp-expm)/(expp+expm)))^0.5;
#
# Return solution
  return(c(ua));
  }