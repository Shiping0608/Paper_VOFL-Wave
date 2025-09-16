function test()
% -------------------------------------------------------------------------
% test()  Driver function for running timesplitting scheme 
%         temporal convergence tests.
%
%   This script sets up a simple test problem and calls the TS solver 
%   for a range of time steps. The purpose is to study convergence 
%   of the scheme for a fixed spatial grid.
%
%   Solving the wave problem:
%
%       u_{tt} = \alpha (-\Delta)^{s(x)} u(x,t) + \lambda f(u),
%
%   with boundary conditions:
%
%       u(x,0) = \phi(x), u_{t}(x,0) = \psi(x).
%
%   Input:
%       None.
%
%   Output:
%       Errors and convergence rates data saved.
%
%   parameters:
%       a       - Half-length of the spatial domain [-a, a]
%       xb, xe  - Left and right boundaries of the domain
%       alpha   - Diffusion coefficient
%       s(x)    - Spatial dependent order
%       T0, T1  - Initial and final time
%       lambda  - Nonlinear term coefficient
%       dts     - Array of time step sizes to test
%       N       - Number of spatial grid points
%
%   Subroutines:
%       splitting(xb,xe,alpha,T0,T1,N,dt,x,s,phi,psi,lambda)
%           – Time stepping solver using Crank–Nicolson scheme.
%
%       error_ref()
%           – Reference error computation routine.
%
%   Author:
%       Shiping Zhou
%
%   Date last updated:
%       September 15, 2025 (Monday)
%
%   References:
%       - Y. Zhang, X. Zhao, and S. Zhou, "Fourier pseudospectral methods 
%         for the variable-order space fractional wave equations"
%         (under review, to be updated)
%
%   Licensing:
%    This code is distributed under the MIT license.
%
% -------------------------------------------------------------------------

a = 32;
xb = -a;
xe = a;
yb = -a;
ye = a;
kappa = 0.2;
T0 = 0;
T1 = 10;

dt = 1e-4;
NN = 2048;
M = 30;
for N = NN
    h = (xe-xb)/N;
    x = [xb:h:xe-h];
    y = [yb:h:ye-h];
    [xx,yy] = meshgrid(x,y);
    s = 1-0.4*cos(0.25*pi*xx).*cos(0.25*pi*yy);
    sp = 0.5*(max(max(s))+min(min(s)));
    
    % initial condition phi and psi
    phi = 3*exp(-5*(xx.^2+yy.^2));
    psi = zeros(size(xx));
    
    %
    lambda = 0; % linear case
    splitting2d(xb,xe,yb,ye,xx,yy,kappa,T0,T1,N,dt,s,sp,phi,psi,lambda,M-1);
end
