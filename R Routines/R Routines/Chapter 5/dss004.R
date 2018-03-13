  dss004=function(xl,xu,n,u) {
#
# An extensive set of documentation comments detailing the derivation
# of the following fourth order finite differences (FDs) is not given
# here to conserve space.  The comments are available in the Matlab 
# version of dss004 in http://www.pdecomp.net.  The derivation is also
# detailed in Schiesser, W. E., The Numerical Method of Lines Integration
# of Partial Differential Equations, Academic Press, San Diego, 1991.
#
# Preallocate arrays
  ux=rep(0,n);
#
# Grid spacing
  dx=(xu-xl)/(n-1);
#
# 1/(12*dx) for subsequent use
  r12dx=1/(12*dx);
#
# ux vector
#
# Boundaries (x=xl,x=xu)
  ux[1]=r12dx*(-25*u[1]+48*u[  2]-36*u[  3]+16*u[  4]-3*u[  5]);
  ux[n]=r12dx*( 25*u[n]-48*u[n-1]+36*u[n-2]-16*u[n-3]+3*u[n-4]);
#
# dx in from boundaries (x=xl+dx,x=xu-dx)
  ux[  2]=r12dx*(-3*u[1]-10*u[  2]+18*u[  3]-6*u[  4]+u[  5]);
  ux[n-1]=r12dx*( 3*u[n]+10*u[n-1]-18*u[n-2]+6*u[n-3]-u[n-4]);
#
# Interior points (x=xl+2*dx,...,x=xu-2*dx)
  for(i in 3:(n-2))ux[i]=r12dx*(-u[i+2]+8*u[i+1]-8*u[i-1]+u[i-2]);
#
# All points concluded (x=xl,...,x=xu)
  return(c(ux));
}
