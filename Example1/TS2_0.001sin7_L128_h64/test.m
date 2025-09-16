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

% -------------------------------------------------------------------------
a = 128;
xb = -a;
xe = a;
alpha = -1;
T0 = 0;
T1 = 40;

dt = 1e-4;
NN = [16384];

M = 25;
for N = NN
    h = (xe-xb)/N;
    x = [xb:h:xe-h]';

    % ---------------------------------------------------------------------
    s = 1+0.001*(1-sin(pi*x/64)).^7;
    sp = 0.5*(max(s)+min(s));

    % ---------------------------------------------------------------------
    % initial conditions
    a = 3; b = 3;
    phi = sech(a*(x+10)-b*T0)+sech(a*(x-10)+b*T0);
    psi = zeros(size(x));
    for i = 1:length(x)
        if abs(x(i)) < 128
            psi(i) = - (b*sinh(b*T0 - a*(x(i) + 10)))./cosh(b*T0 - a*(x(i) + 10)).^2 ...
                     - (b*sinh(b*T0 + a*(x(i) - 10)))./cosh(b*T0 + a*(x(i) - 10)).^2;
        end
    end
    splitting(xb, xe, N, x, T0, T1, dt, alpha, phi, psi, sp, s, M-1)
end
