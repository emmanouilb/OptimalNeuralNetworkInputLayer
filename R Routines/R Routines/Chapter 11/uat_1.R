  uat_1=function(x,t){
#
# Function uat_1 computes the time derivative of the 
# analytical solution of the Boussinesq equation
#
# Analytical solution derivative
  xi=(1-c^2)^(0.5)/2*(x-c*t);
  uat=(3/2)*c*(1-c^2)^(3/2)*sinh(xi)*(cosh(xi))^(-3);
  return(c(uat));
  }