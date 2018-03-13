function [D0,D2,D4,Mt,Dl,Dr]=init(nx,dx,dt,mu,alpha)
% this  function would only be needed for Crank-Nicholsn which requires a number
% of matrices that need to be generated initially (and whenever dx etc. gets changed)
% for forward Euler init.m is not needed.

e=ones(nx,1);

D0=spdiags([e],0,nx,nx);

