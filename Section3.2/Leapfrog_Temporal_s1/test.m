function test()
% -------------------------------------------------------------------------
% test()  Driver function for running Leapfrog scheme 
%         temporal convergence tests.
%
%   This script sets up a simple test problem and calls the LP solver 
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
%       Leapfrog(xb,xe,alpha,T0,T1,N,dt,x,s,phi,psi,lambda)
%           – Time stepping solver using Crank–Nicolson scheme.
%
%       error_ref()
%           – Reference error computation routine.
%
%   Author:
%       Shiping Zhou
%
%   Date last updated:
%       [2025-09-15]
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

% -------------------------------------------------------------------------
% Setup parameters
a = 32;
xb = -a;
xe = a;
alpha = -1;
T0 = 0;
T1 = 1;

dts = [1/128 1/256 1/512 1/1024 1/2048 1e-4];
N = 4096;
M = 40;
for dt = dts
    h = (xe-xb)/N;
    x = [xb:h:xe-h]';
    
    % ---------------------------------------------------------------------
    s = 1+zeros(size(x));
    
    % ---------------------------------------------------------------------
    % initial conditions
    phi = exp(-x.^2);
    psi = zeros(size(x));
    
    % ---------------------------------------------------------------------
    % nonlinear term on the right-hand side lambda*u^3
    lambda = 1;
    
    Leapfrog(xb,xe,alpha,T0,T1,N,dt,x,s,phi,psi,lambda,M);
end
error_ref()
