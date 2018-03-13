function [U,x,t] = ksintegrate(u, Lx, dt, Nt, nplot);
% KS eqn: u_t = -u*u_x - u_xx - u_xxxx, periodic BCs 


Nx = length(u);                % number of gridpoints
kx = [0:Nx/2-1 0 -Nx/2+1:-1];  % integer wavenumbers: exp(2*pi*i*kx*x/L)
alpha = 2*pi*kx/Lx;            % real wavenumbers:    exp(i*alpha*x)
D = i*alpha;                   % D = d/dx operator in Fourier space
L = alpha.^2 - alpha.^4;       % linear operator -D^2 - D^3 in Fourier space
G = -0.5*D;                    % -1/2 D operator in Fourier space
NT = floor(Nt/nplot) + 1;      % number of saved time steps

x = (0:Nx-1)*Lx/Nx;
t = (0:NT)*dt*nplot;

% some convenience variables
dt2  = dt/2;
dt32 = 3*dt/2;
A =  ones(1,Nx) + dt2*L;
B = (ones(1,Nx) - dt2*L).^(-1);

Nn  = G.*fft(u.*u); % -u u_x, spectral
Nn1 = Nn;        

U(1,:) = u;   % save u(x,0), physical
u = fft(u);   % u, spectral
np = 2;       % n saved,

% timestepping loop
for n = 1:Nt

  Nn1 = Nn;   % shift nonlinear term in time: N^{n-1} <- N^n
  Nn  = G.*fft(real(ifft(u)).^2); % compute Nn = -u u_x

  u = B .* (A .* u + dt32*Nn - dt2*Nn1);

  if (mod(n, nplot) == 0)
    U(np, :) = fft(u);
    np = np + 1;
  end

end