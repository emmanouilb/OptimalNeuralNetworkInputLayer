  u2a=function(x,t){
#
# Function u2a computes the analytical solution 
# for u2(x,t) of the Keller-Segel equations
#
  z=x-c*t;
  u2a=(c^2/(k*D))*exp(-c*z/D)/(1+exp(-c*z/D))^2;
#
# Return solution
  return(c(u2a));
  }