function CN(xb, xe, alpha, T0, T1, N, x, s, sp, dt, lambda, M)
% -------------------------------------------------------------------------
% CN() Solving the variable-order fractional wave equation with
%      Crank-Nicolson for temporal discretization
%
%   Input:
%       xb, xe   - Left and right spatial boundaries
%       alpha    - Diffusion coefficient
%       T0, T1   - Initial and final time
%       N        - Number of spatial grid points
%       x        - Spatial grid points (column vector)
%       s, sp    - Vector for s(x) and its average value
%       dt       - Time step size
%       lambda   - Nonlinear term coefficient
%       M        - Number of modes / truncation parameter for LMu
%
%   Output:
%       Solution saved for computing errors.
%
%   Subroutines:
%       LMu(u_hat, M, mu, s, sp) 
%           – Evaluates the fractional operator in Fourier space.
%
%       CG_Solver(alpha, dt, u_n, f, M, mu, s, sp, lambda) 
%           – Solves the implicit CN system using conjugate gradient method.
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



NT = round((T1-T0)/dt);
mu = [0:N/2-1,-N/2:-1]*(2*pi/(xe-xb));

% solve the nonlinear system
u0 = exp(-x.^2);     % u(x,0)
v0 = zeros(size(x)); % u_{t}(x,0)

u0_hat = fft(u0);
LMu0 = LMu(u0_hat, M, mu', s, sp);
u1 = u0 + dt*v0 + 0.5*dt^2*(alpha*LMu0+lambda*u0.^3);

u_nm1_hat = u0_hat;
u_nm1 = u0;
u_n = u1;
for n = 1:NT-1
    tnp1 = (n+1)*dt;
    LMu_nm1 = LMu(u_nm1_hat, M, mu', s, sp);  %D*u^{n-1}
    f = 0.5*dt^2*(alpha*LMu_nm1+lambda*u_nm1.^3) - u_nm1 + 2*u_n;
    [unp1, ~, ~] = CG_Solver(alpha, dt, u_n, f, M, mu', s, sp, lambda);
    u_nm1 = u_n;
    u_nm1_hat = fft(u_nm1);
    u_n = unp1;
end

save("sol_sin_tau"+num2str(NT)+"_t"+num2str(double(tnp1)) ...
             +".mat",'tnp1','x','unp1','N');