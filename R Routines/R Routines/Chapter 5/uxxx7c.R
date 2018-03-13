  uxxx7c=function(xl,xu,n,u){
#
# Function uxxx7c computes the derivative uxxx 
# in the KdV equation
#
# Coefficient
  r8dx3=1/(8*(dx^3));
#
# uxxx
  uxxx=rep(0,n);
  for(i in 1:n){
#
#   At the left end, uxxx = 0
    if(i<4){uxxx[i]=0;} 
#
#   At the right end, uxxx = 0
    if(i>(n-3)){uxxx[i]=0;} 
#
#   Interior points
    if((i>=4)&(i<=(n-3))){
    uxxx[i]=r8dx3*
       (   1*u[i-3]- 
           8*u[i-2]+
          13*u[i-1]+
           0*u[i  ]-
          13*u[i+1]+
           8*u[i+2]-
           1*u[i+3]);
    }
  }
#
# Return derivative
  return(c(uxxx));
  } 