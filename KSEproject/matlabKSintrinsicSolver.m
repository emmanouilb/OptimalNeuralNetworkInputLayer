  % from http://chaosbook.org/extras/KSEproject/html/index.html

  N = 22;  L = 22;  h = 0.25;  nstp = 3000;
  a0 = zeros(N-2,1);  a0(round(size(a0)/4*1):round(size(a0)/4*3)) = 0.00001; % just some initial condition
  [tt, at] = ksfmstp(a0, L, h, nstp, 1);
  %fig1 = figure('pos',[0 0 600 400],'color','w'); plot(tt,at,'.-');
  %title(strcat('Solution of Kuramoto-Sivashinsky equation with L = ',int2str(N),': Fourier modes'));
  %xlabel('Time'); ylabel('Fourier modes');
  
  [x, ut] = ksfm2real(at, L);
  fig2 = figure('pos',[5 270 600 200],'color','w');
  pcolor(tt,x,ut); shading interp; caxis([-3 3]);
  title(strcat('Solution u(x,t) of Kuramoto-Sivashinsky equation, system size L = ' ,num2str(L)));
  xlabel('Time'); ylabel('x','rotat',0);
  display(size(ut));
  csvwrite('outputKSL=22.csv',ut)
  csvwrite('outputKSL=22Square.csv',ut.^2)
  csvwrite('outputKSL=22Absolutee.csv',abs(ut))
