function fractionalWave_LFFS(xb, xe, yb, ye, xx, yy, T0, T1, N, dt, ...
                             s1, sp1, s2, sp2, phi, psi, c, eta, tau, M)
% =========================================================================
% Function: fractionalWave_LFFS
%
% Purpose:
%   Solve the two-dimensional fractional wave equation using a leapfrog
%   time discretization scheme with a matrix-free acceleration
%   approach.
%
% Syntax:
%   fractionalWave_LFFS(xb, xe, yb, ye, xx, yy, T0, T1, N, dt, ...
%                       s1, sp1, s2, sp2, phi, psi, c, eta, tau)
%
% Input:
%   - xb, xe   : domain boundaries in x-direction (scalars)
%   - yb, ye   : domain boundaries in y-direction (scalars)
%   - xx, yy   : 2D spatial meshgrid arrays for x and y coordinates
%   - T0, T1   : start and end times (scalars)
%   - N        : number of spatial grid points in each direction (integer)
%   - dt       : time step size (scalar)
%   - s1, sp1  : variable-order and its average for the first term
%   - s2, sp2  : variable-order and its average for the second term
%   - phi      : initial displacement field, u(x,y,0)
%   - psi      : initial velocity field, u_t(x,y,0)
%   - c        : velocity profile c(x,y)
%   - eta      : coefficient function η(x,y)
%   - tau      : coefficient function τ(x,y)
%
% Output:
%   - Saves solution at final time tn into a MAT file:
%       "data_RickerWavelet_N<grid>_t<time>.mat"
%     containing variables {xx, yy, tn, u}.
%
% Description:
%   This function computes the solution to the 2D variable-order fractional
%   wave equation:
%
%       (1/c(x,y)^2) * u_tt =
%           η(x,y)(-Δ)^{1+s(x,y)} u(x,y,t)
%         + τ(x,y)(-Δ)^{1/2+s(x,y)} u(x,y,t),
%
%   with initial displacement φ(x,y) and initial velocity ψ(x,y).
%   The scheme uses:
%     - Leapfrog time discretization.
%     - Matrix-free acceleration for evaluating fractional Laplacians.
%     - Fast spectral approximation (M terms) for efficient evaluation.
%
% Example:
%   % Define grid and parameters externally, then call:
%   fractionalWave_LFFS(0,2,0,2,xx,yy,0,0.5,128,1e-4,...
%                       s1,sp1,s2,sp2,phi,psi,c,eta,tau);
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

NT = round((T1-T0)/dt);


% differential matrix D2
mu_l = [0:N/2-1,-N/2:-1]*(2*pi/(xe-xb));
mu_m = [0:N/2-1,-N/2:-1]*(2*pi/(ye-yb));
[mu_ll, mu_mm] = meshgrid(mu_l,mu_m);
mu_lm = sqrt(mu_ll.^2+mu_mm.^2);

u0 = phi;
v0 = psi;
 
u0_hat = fft2(u0);


% get u1 by 1st-order Taylor expansion
DM1s_u0 = LMu(u0_hat, M, mu_lm, s1, sp1);
u1 =  u0 + dt*v0 + 0.5*dt^2*c.^2.*eta.*DM1s_u0;
u1_hat = fft2(u1);

DM2s_unm2 = LMu(u0_hat, M, mu_lm, s2, sp2);
DM2s_unm1 = LMu(u1_hat, M, mu_lm, s2, sp2);


unm1 = u1; % u_{n-1}
unm2 = u0; % u_{n-2}
unm1_hat = fft2(double(unm1));
for n = 2:NT
    tn = n*dt;
    % ---------------------------------------------------------------------
    % step one, half timestep, 0.5*dt
    DM1s_unm1 = LMu(unm1_hat, M, mu_lm, s1, sp1);
    un = 2*unm1 + dt^2*c.^2.*eta.*DM1s_unm1 ...
                + dt*c.^2.*tau.*(DM2s_unm1-DM2s_unm2) - unm2;
    %
    unm2 = unm1;
    unm1 = un;
    DM2s_unm2 = DM2s_unm1;
    unm1_hat = fft2(un);
    DM2s_unm1 = LMu(unm1_hat, M, mu_lm, s2, sp2);
end
u = un;
save("data_RickerWavelet_N"+num2str(N)+"_t"+num2str(tn)+".mat", ...
             'xx','yy','tn','u');
end