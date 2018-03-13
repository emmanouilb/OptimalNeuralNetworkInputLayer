  ua=function(x,t){
#
# Funtion ua computes the exact solution of the KdV 
# equation for comparison with the numerical solution. 
#
# Analytical solution
  expm=exp(-1/2*sqrt(c1)*(x-c1*t));
  expp=exp( 1/2*sqrt(c1)*(x-c1*t)); 
  ua=(1/2)*c1*(2/(expp+expm))^2;
#
# Return analytical solution
  return(c(ua));
}