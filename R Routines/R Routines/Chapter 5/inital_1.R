  inital_1=function(x,t0){
#
# Function inital_1 sets the initial condition for the 
# KdV equation
#
  u0=rep(0,n);
#
# Case 1 - single pulse
  if(ncase==1){
    for(i in 1:n){
      u0[i]=ua(x[i],0);
    }
  } 
#
# Case 2 - two pulses
  if(ncase==2){
    for(i in 1:n){
      expm=exp(-1/2*sqrt(c1)*(x[i]+15));
      expp=exp( 1/2*sqrt(c1)*(x[i]+15)); 
      pulse1=(1/2)*c1*(2/(expp+expm))^2;  
      expm=exp(-1/2*sqrt(c2)*(x[i]-15));
      expp=exp( 1/2*sqrt(c2)*(x[i]-15)); 
      pulse2=(1/2)*c2*(2/(expp+expm))^2;      
      u0[i]=pulse1+pulse2;
    }
  }
#
# Return IC
  return(c(u0));
  }