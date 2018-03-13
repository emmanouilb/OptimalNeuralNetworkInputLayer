  uat_1=function(x,t){
#
# Function uat_1 computes the time derivative of the 
# exact solution of the cubic Klein-Gordon equation
#
# Analytical solution derivative
  uat=B*(1/cos(K*(x+c*t))^2*(K*c));
  return(c(uat));
  }