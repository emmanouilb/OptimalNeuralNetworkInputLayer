  simp=function(xl,xu,n,u){
#
# Function simp computes three integral invariants 
# by Simpson's rule
#
  uint=rep(0,3);
  h=(xu-xl)/(n-1); 
#
  for(int in 1:3){
#
#   Conservation of mass
    if(int==1){  
       uint[1]=u[1]-u[n]; 
       iv=seq(from=3,to=n,by=2);
       for(i in iv){
         uint[1]=uint[1]+4*u[i-1]+2*u[i];  
       }
       uint[1]=h/3*uint[1];
    } 
#
#   Conservation of energy
    if(int==2){  
       uint[2]=u[1]^2-u[n]^2; 
       iv=seq(from=3,to=n,by=2);
       for(i in iv){
         uint[2]=uint[2]+4*u[i-1]^2+2*u[i]^2;  
       }
       uint[2]=(1/2)*h/3*uint[2];
    }   
#
#   Whitham conservation
    if(int==3){
       ux=dss004(xl,xu,n,u); 
       uint[3]=2*u[1]^3-ux[1]^2-(2*u[n]^3-ux[n]^2); 
       iv=seq(from=3,to=n,by=2);
       for(i in iv){
         uint[3]=uint[3]+4*(2*u[i-1]^3-ux[i-1]^2)+
                         2*(2*u[i]^3  -ux[i]^2);  
       }
       uint[3]=h/3*uint[3];
    }       
  } 
#
# Return vector of integrals
  return(c(uint));
  } 