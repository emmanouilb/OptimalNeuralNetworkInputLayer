  u1a=function(x,t){
#
# Function u1a computes the analytical solution 
# for u1(x,t) of the Keller-Segel equations
  z=x-c*t;
  u1a=1/(1+exp(-c*z/D));
#
# Return solution
  return(c(u1a));
  }