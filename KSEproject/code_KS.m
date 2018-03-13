% Newton-CG method for time-space-periodic solutions
% in the KS equation: u_t+uu_x+u_{xx}+gamma*u_{xxxx}=0.
% In this code, z represents scaled time tau in the paper.

gamma=1;
Nx=200; Nz=200; Lx=2*pi; Lz=2*pi; errormax=1e-8; errorCG=1e-4;
dx=Lx/Nx; x=0:dx:Lx-dx;  kx=[0:Nx/2-1  -Nx/2:-1]*2*pi/Lx;
dz=Lz/Nz; z=0:dz:Lz-dz;  kz=[0:Nz/2-1  -Nz/2:-1]*2*pi/Lz;
[X,Z]=meshgrid(x,z); [KX,KZ]=meshgrid(kx,kz);
KX2=-KX.*KX+gamma*KX.^4;
u0=-7*sin(3*X)-3*sin(Z).*(sin(4*X)-sin(5*X))-cos(Z).*sin(X);
u=u0;    % i.c.

tic;
nnt=0;   % nnt: # of Newton steps
ncg=0;   % ncg: # of CG iterations
while 1                        % Newton-CG iterations start
  nnt=nnt+1;
  ufft=fft2(u);
  F=-real(u.*ifft2(i*KX.*ufft)+ifft2(KX2.*ufft));
  uz=real(ifft2(i*KZ.*ufft));
  omega=sum(sum(uz.*F))/sum(sum(uz.*uz));
  L0u=omega*uz-F;
  uerror(nnt)=max(max(abs(L0u))); uerror(nnt)
  numcg(nnt)=ncg; time(nnt)= toc;
  if uerror(nnt) < errormax
    break
  end

  P=@(W)  real(ifft2(( omega*i*KZ+KX2).*fft2(W)) ...
               +ifft2(i*KX.*fft2(u.*W)));
  PA=@(W) real(ifft2((-omega*i*KZ+KX2).*fft2(W)) ...
               -u.*ifft2(i*KX.*fft2(W)));

  c=30; fftM=omega^2*KZ.*KZ+KX2.*KX2+c;    % Preconditioner
  du=0*Z;                             % CG iterations start
  R=-PA(L0u);
  MinvR=real(ifft2(fft2(R)./fftM));
  R2=sum(sum(R.*MinvR)); R20=R2;
  D=MinvR;
  while (R2 > R20*errorCG^2)
     PD=P(D);
     L1D=PD-sum(sum(uz.*PD))/sum(sum(uz.*uz))*uz;
     PAL1D=PA(L1D);
     a=R2/sum(sum(D.*PAL1D));
     du=du+a*D;
     R=R-a*PAL1D;
     MinvR=real(ifft2(fft2(R)./fftM));
     R2old=R2;
     R2=sum(sum(R.*MinvR));
     b=R2/R2old;
     D=MinvR+b*D;
     ncg=ncg+1;
  end                                   % CG iterations end
  u=u+du;
end                              % Newton-CG iterations end

% plotting of numerical results
subplot(221); imagesc(x, [z z+Lz], [u0; u0]);
axis xy; colorbar;
xlabel('x'); ylabel('\tau','rotation',0); title('(a)');
subplot(222); imagesc(x, [z z+Lz]/omega, [u; u]);
axis xy; colorbar;
xlabel('x'); ylabel('t','rotation',0); title('(b)');
subplot(223); semilogy(numcg, uerror, numcg, uerror, 'o');
xlabel('number of CG iterations'); ylabel('solution error');
title('(c)');
subplot(224); semilogy(time, uerror, time, uerror, 'o');
xlabel('time (seconds)'); ylabel('solution error');
title('(d)');
format long; period=2*pi/omega