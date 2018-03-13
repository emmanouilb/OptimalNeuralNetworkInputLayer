% KS with Euler periodic boundary conditions

function ks

[nx,xmax,tmax,dt,dx,D]=default;

ind=0;
while ind < 98

ind=input(' 0: compute 1: i.c. 2: nx,xmax,tmax,dt 3: D 99: quit ');

 if ind == 0

   [y]=go(x,y,t,tmax,xmax,nx,dt,D);

 end

 if ind == 1

% for Crank-Nicholson init determines the various matrices
%   [D0,D2,D4,Mt,Dl,Dr]=init(nx,dx,dt,mu,alpha);
   [t,x,y]=ic(nx,dx,xmax);

 end

 if ind == 2

  nx
  xmax
  tmax
  dt
  nx=input(' nx =');
  xmax=input(' xmax =');
  tmax=input(' tmax =');
  dt=input(' dt =');
  dx=xmax/nx;
%  [D0,D2,D4,Mt,Dl,Dr]=init(nx,dx,dt,mu,alpha);


 end

 if ind == 3

  D
  D=input(' D =');
%  [D0,D2,D4,Mt,Dl,Dr]=init(nx,dx,dt,mu,alpha);

 end

end
 
