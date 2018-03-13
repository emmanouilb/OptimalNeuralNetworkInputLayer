function [y]=go(x,y,t,tmax,xmax,nx,dt,D)

nt=0;
ntplot=10;
fig1=figure(1);
lin1=line(x,y,'EraseMode','background');
axis([0 xmax -1 1]);
drawnow;

dx=xmax/nx;

while t < tmax
ym=[y(nx) y(1:nx-1)];
ymm=[ym(nx) ym(1:nx-1)];
yp=[y(2:nx) y(1)];
ypp=[yp(2:nx) yp(1)];
yn=y+dt*(D*(yp-2*y+ym)/dx^2-1/dx^4*(ypp-4*yp+6*y-4*ym+ymm)+...
1/(2*dx)*(yp-ym).*y);
y=yn;

if rem(nt,ntplot) == 0
delete(lin1);
lin1=line(x,y,'EraseMode','background');
drawnow;
end

nt=nt+1;
t=t+dt;
end

