clear;
clc;
Lx = 200;
Nx = 1024;
dt = 1/16;
nplot = 8;
Nt = 1600;

x = Lx*(0:Nx-1)/Nx;
u = abs(cos(x) + 0.1*cos(x/16).*(1+2*sin(x/16)));

[U,x,t] = ksintegrate(u, Lx, dt, Nt, nplot);

size(U)
size(x)
size(t)

realUpart = transpose(real(U));

  fig2 = figure('pos',[5 270 600 200],'color','w');
  pcolor(t(1:size(U,1)),transpose(x),realUpart); shading interp; colormap jet;
  view([-90 90]);
  %caxis([-3 3]);
  title(strcat('Solution u(x,t) of Kuramoto-Sivashinsky equation, system size L = ' ,num2str(Lx)));
  xlabel('Time'); ylabel('x','rotat',0);

  disp(max(max(realUpart)));
  
  disp(min(min(realUpart)));
