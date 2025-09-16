function Leapfrog(xb,xe,alpha,T0,T1,N,dt,x,s,phi,psi,lambda,M)
% -------------------------------------------------------------------------
% Leapfrog() Solving the variable-order fractional wave equation with
%            Leapfrog for temporal discretization
%
%   Input:
%       xb, xe   - Left and right spatial boundaries
%       alpha    - Diffusion coefficient
%       T0, T1   - Initial and final time
%       N        - Number of spatial grid points
%       dt       - Time step size
%       x        - Spatial grid points (column vector)
%       s, sp    - Vector for s(x) and its average value
%       phi, psi - Initial condition u(x,0) and u_{t}(x,0)
%       lambda   - Nonlinear term coefficient
%       M        - Number of modes / truncation parameter for LMu
%
%   Output:
%       Solution saved for computing errors.
%
%   Subroutines:
%       LMu(u_hat, M, mu, s, sp) 
%           – Evaluates the fractional operator in Fourier space with 
%             fast algorithm.
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


% u(x,0)
u0 = phi;

u0_hat = fft(double(u0));
% matrix free version -----------------------------------------------------
sp = ones(length(s),1)*0.5*(max(s)+min(s));

if all(s == s(1))
    LMu0 = ifft(abs(mu').^(2*s).*u0_hat);
else
    LMu0 = LMu(u0_hat, M, mu', s, sp);
end

% u_t(x,0)
v0 = psi;

% get u1 by taylor expansion
u1 = u0 + dt*v0 + 0.5*dt^2*(alpha*LMu0+lambda*u0.^3);


unm1 = u1; % unm1 : u_{n-1}
unm2 = u0; % unm2 : u_{n-2}
for n = 2:NT
    nt = n*dt;
    unm1_hat = fft(double(unm1));
    if all(s == s(1))
        LMunm1 = ifft(abs(mu').^(2*s).*unm1_hat);
    else
        LMunm1 = LMu(unm1_hat, M, mu', s, sp);
    end
    un = 2*unm1 + dt^2*(alpha*LMunm1+lambda*unm1.^3) - unm2;
    unm2 = unm1;
    unm1 = un;
end
save("sol_s0_tau"+num2str(NT)+"_t"+num2str(nt)+".mat",'nt','x','un','N');
