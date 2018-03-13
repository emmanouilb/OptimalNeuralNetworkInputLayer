% Parameters:
m = 1024; n = m;
% number of grid points
h = 1e-1; T = 100;
% time-step and final time
u0 = @(x,y,z) cos(40*x)+cos(40*y)+cos(40*z);
th = pi/8; c = cos(th); s = sin(th);
u0 = 1/3*spherefun(@(x,y,z) u0(c*x-s*z,y,s*x+c*z)); % initial condition
v0 = reshape(coeffs2(u0, m, n), m*n, 1);            % Fourier coefficients
% Nonlinear operator (evaluated in value space):
g = @(u) u - (1+1.5i)*u.*(abs(u).^2);               % N(u) = u-(1+1.5)i*u*|u|^2
c2v = @(u) trigtech.coeffs2vals(u);                 % coeffs to values in 1D
c2v = @(u) c2v(c2v(reshape(u,m,n)).').';            % coeffs to values in 2D
v2c = @(u) trigtech.vals2coeffs(u);                 % values to coeffs in 1D
v2c = @(u) reshape(v2c(v2c(u).').',m*n,1);          % values to coeffs in 2D
N = @(u) v2c(g(c2v(u)));                            % nonlinear operator
% Construct the Laplacian matrix (multiplied by Tsin2 and 1e-4):
Dm = spdiags(1i*[0,-m/2+1:m/2-1]', 0, m, m);
D2m = spdiags(-(-m/2:m/2-1).^2', 0, m, m);
D2n = spdiags(-(-n/2:n/2-1).^2', 0, n, n);
Im = speye(m); In = speye(n);
P = speye(m+1); P = P(:, 1:m); P(1,1) = .5; P(m+1,1) = .5;
Q = speye(m+1+4); Q = Q(3:m+2,:); Q(1,3) = 1; Q(1,m+3) = 1;
Msin2 = toeplitz([1/2, 0, -1/4, zeros(1, m+2)]);
Msin2 = sparse(Msin2(:, 3:m+3));
Tsin2 = round(Q*Msin2*P, 15);
Mcossin = toeplitz([0, 0, 1i/4, zeros(1, m+2)]);
Mcossin = sparse(Mcossin(:, 3:m+3));
Tcossin = round(Q*Mcossin*P, 15);
Lap = 1e-4*(kron(In, Tsin2*D2m + Tcossin*Dm) + kron(D2n, Im));
% Compute LU factorizations of LIRK4 matrices:
Tsin2 = kron(In, Tsin2);
[L, U] = lu(Tsin2); [La, Ua] = lu(Tsin2 - 1/4*h*Lap);
% Time-stepping loop:
itermax = round(T/h); v = v0;
for iter = 1:itermax
Nv = N(v); w = Tsin2*v;
wa = w + h*Tsin2*1/4*Nv;
a = Ua\(La\wa); Na = N(a);
wb = w + h*Lap*1/2*a + h*Tsin2*(-1/4*Nv + Na);
b = Ua\(La\wb); Nb = N(b);
wc = w + h*Lap*(17/50*a - 1/25*b) + h*Tsin2*(-13/100*Nv + 43/75*Na + 8/75*Nb);
c = Ua\(La\wc); Nc = N(c);
wd = w + h*Lap*(371/1360*a - 137/2720*b + 15/544*c) ...
+ h*Tsin2*(-6/85*Nv + 42/85*Na + 179/1360*Nb - 15/272*Nc);
d = Ua\(La\wd); Nd = N(d);
we = w + h*Lap*(25/24*a - 49/48*b + 125/16*c - 85/12*d) ...
+ h*Tsin2*(79/24*Na - 5/8*Nb + 25/2*Nc - 85/6*Nd);
e = Ua\(La\we); Ne = N(e);
v = v + h*(U\(L\(Lap*(25/24*a - 49/48*b + 125/16*c - 85/12*d +1/4*e)))) ...
+ h*(25/24*Na - 49/48*Nb + 125/16*Nc - 85/12*Nd + 1/4*Ne);
end
vals = c2v(v);
% tramsform to value space
vals = vals([m/2+1:m 1], :);                        % restrict to [-pi,pi]x[0,pi]
u = spherefun(real(vals)); plot(u)                  % output real(u) an