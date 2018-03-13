function [t,x,y]=ic(nx,dx,xmax)

nwave=input(' nwave: (0=random)');

x=[1:nx]*dx;

y=.5*sin(2*pi*nwave*x/xmax);
   
if nwave == 0
 % random initial conditions
 amp=1
 y=amp*randn(1,nx);
end

t=0;

clf;