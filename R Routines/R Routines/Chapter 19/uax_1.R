  uax_1=function(x,t){
#
# Function uax_1 computes the derivative of two 
# analytical solutions of the Kuramoto-Sivashinsky 
# equation
#
# Derivative of analytical solution
  if(ncase==1){
    F=k/(1+exp(-k*x-lambda*t));
    Fx=(k^2)*(1+exp(-k*x-lambda*t))^(-2)*exp(-k*x-lambda*t);
    uax=c1*Fx+2*c2*F*Fx+3*c3*F^2*Fx;
  }
  if(ncase==2){
    arg=0.5*k*x-(15/19)*k^2*t;
    H=tanh(arg);
    Hx=(1-tanh(arg)^2)*0.5*k;
    uax=(15/19)*k*(33*H^2*Hx-9*Hx);
  } 
#
# Return analytical solution derivative            
  return(c(uax))
}