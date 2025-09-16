function test()
% TEST   Driver function for 2D fractional wave equation 
%        in a heterogeneous medium.
%
%   Consider a two-layer heterogeneous medium on the computational domain
%       Ω = (0, 2)^2.
%
%   Governing equation:
%       (1 / c(x,y)^2) * u_tt(x,y,t) =
%           η(x,y) * (-Δ)^(1+s(x,y)) u(x,y,t) +
%           ζ(x,y) * (-Δ)^(1/2+s(x,y)) u(x,y,t),
%
%   where
%       s(x,y) = 1 + a₁ + a₂ * tanh(100(y - 1)),
%
%   Initial condition:
%       u(x,y,0) = (1 - 2π² f² ((x - x_c)² + (y - y_c)²)) ...
%                   * exp(-π² f² ((x - x_c)² + (y - y_c)²)),
%       with f = 25 and (x_c, y_c) = (1, 0.85).
%
%   Coefficients:
%       η(x,y)   = -c(x,y)^(-2s - 2) * cos(πs),
%       ζ(x,y)   = -c(x,y)^(-2s - 1) * sin(πs).
%
%   Input:
%      None;
%
%   Output:
%       Solution saved.
%
%   Author:
%       Shiping Zhou
%
%   Date latest update:
%       December 30, 2023 (Saturday)
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

L = 2;
disp("Computational domain is [0,"+num2str(L)+"]^2")
xb = 0;
xe = L;
yb = 0;
ye = L;
T0 = 0;
T1 = 0.5;

dt = 5e-5;
NN = 1000;
M = 20;

for N = NN
    h = (xe-xb)/N;
    x = [xb:h:xe-h];
    y = [yb:h:ye-h];
    [xx,yy] = meshgrid(x,y);
    interface = L/2;
    a1 = 0.01;
    a2 = 0.02;
    s = a1+a2/2 + a2/2*tanh(100*(yy-interface));
    s1 = 1+s;
    sp1 = 1+0.5*(max(max(s))+min(min(s)));
    s2 = 0.5+s;
    sp2 = 0.5+0.5*(max(max(s))+min(min(s)));
    % initial condition phi and psi
    f = 25;
    x_c = L/2;
    y_c = interface*0.85;
    phi = (1-2*pi^2*f^2*((xx-x_c).^2+(yy-y_c).^2)) ...
               .*exp(-pi^2*f^2*((xx-x_c).^2+(yy-y_c).^2));
    psi = zeros(size(xx));
    
    % initial velocity
    v0_up = 2.2/3.6;
    v0_bottom = 1;
    c0 = v0_up*(yy>=interface)+v0_bottom*(yy<interface);
    c = cos(0.5*pi*s); % velocity
    %w0 = f; % reference requency, same as the source
    eta = -c0.^(-2*s-2).*cos(pi*s);
    tau = -c0.^(-2*s-1).*sin(pi*s);
    %
    fractionalWave_LFFS(xb,xe,yb,ye,xx,yy,T0,T1,N,dt,s1,sp1,s2,sp2, ...
                        phi,psi,c,eta,tau,M);
end
