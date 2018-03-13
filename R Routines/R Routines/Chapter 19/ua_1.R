  ua_1=function(x,t){
#
# Function ua_1 computes two analytical solutions
# of the Kuramoto-Sivashinsky equation
#
# Analytical solution
  if(ncase==1){
    F=k/(1+exp(-k*x-lambda*t));
    ua=c0+c1*F+c2*F^2+c3*F^3;
  }
  if(ncase==2){
    H=tanh(0.5*k*x-(15/19)*k^2*t);
    ua=(15/19)*k*(11*H^3-9*H+2);
  } 
#
# Return analytical solution            
  return(c(ua))
}