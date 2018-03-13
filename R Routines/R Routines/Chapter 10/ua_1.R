  ua_1=function(x,t){
#
# Function ua_1 computes the exact solution of the cubic
# Klein-Gordon equation for comparison with the numerical 
# solution 
#
# Analytical solution
  ua=B*tan(K*(x+c*t));
  return(c(ua));
  }