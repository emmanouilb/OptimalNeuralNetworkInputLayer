  dss044=function(xl,xu,n,u,ux,nl,nu) {
#
# An extensive set of documentation comments detailing the derivation
# of the following fourth order finite differences (FDs) is not given
# here to conserve space.  The comments are available in the Matlab 
# version of dss044 in http://www.pdecomp.net.  The derivation is also
# detailed in Schiesser, W. E., The Numerical Method of Lines Integration
# of Partial Differential Equations, Academic Press, San Diego, 1991.
#
# Preallocate arrays
  uxx=rep(0,n);
#
# Grid spacing
  dx=(xu-xl)/(n-1);
#
# 1/(12*dx**2) for subsequent use
  r12dxs=1/(12*dx^2);
#
# uxx vector
#
# Boundaries (x=xl,x=xu)
  if(nl==1)
    uxx[1]=r12dxs*
           (45*u[  1]-154*u[  2]+214*u[  3]-
           156*u[  4] +61*u[  5] -10*u[  6]);
  if(nu==1)
    uxx[n]=r12dxs*
           (45*u[  n]-154*u[n-1]+214*u[n-2]-
           156*u[n-3] +61*u[n-4] -10*u[n-5]);
  if(nl==2)
    uxx[1]=r12dxs*
           (-415/6*u[  1] +96*u[  2]-36*u[  3]+
              32/3*u[  4]-3/2*u[  5]-50*ux[1]*dx);
  if(nu==2)
    uxx[n]=r12dxs*
           (-415/6*u[  n] +96*u[n-1]-36*u[n-2]+
              32/3*u[n-3]-3/2*u[n-4]+50*ux[n]*dx);
#
# dx in from boundaries (x=xl+dx,x=xu-dx)
    uxx[  2]=r12dxs*
           (10*u[  1]-15*u[  2]-4*u[  3]+14*u[  4]-
             6*u[  5]   +u[  6]);
    uxx[n-1]=r12dxs*
           (10*u[  n]-15*u[n-1]-4*u[n-2]+14*u[n-3]-
             6*u[n-4]   +u[n-5]);
#
# Remaining interior points (x=xl+2*dx,...,x=xu-2*dx)
  for(i in 3:(n-2))
    uxx[i]=r12dxs*
           (-u[i-2]+16*u[i-1]-30*u[i]+16*u[i+1]-u[i+2]);
#
# All points concluded (x=xl,...,x=xu)
  return(c(uxx));
}
