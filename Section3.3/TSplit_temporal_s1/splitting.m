function splitting(xb, xe, N, T0, T1, dt, alpha, phi, psi, sp, s, M)
% -------------------------------------------------------------------------
% splitting() Solving the variable-order fractional wave equation with
%            timesplitting for temporal discretization
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



NT = round((T1-T0)/dt); % number of time levels
h = (xe-xb)/N;
x = [xb:h:xe-h]; % mesh

mu = [0:N/2,-N/2+1:-1]*(2*pi/(xe-xb));
% DM2_sp = real((exp(1i*x.'*mu).*(abs(mu).^(2*sp+2)))*exp(-1i*mu.'*x)/N);
% DM2_s = real((exp(1i*x.'*mu).*(abs(mu).^(2*s+2)))*exp(-1i*mu.'*x)/N);

wl = sqrt(-alpha)*abs(mu').^(sp(1)+1);
emu = exp(1i*(x-xb)'*mu);

u0 = phi; % u(x,0)
v0 = psi; % u_{t}(x,0)

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% need to write fft and ifft with higher precision, maybe vpa with DFT
u0_hat = fft(u0);
v0_hat = fft(v0);

s0 = 0.5*(max(s)+min(s));

c1 = 0.5;
c2 = 1;
lambda = 1;

for n=1:NT
    % ---------------------------------------------------------------------
    % step one, half timestep, 0.5*dt
    temp1 = u0_hat.*cos(c1*dt*wl)+sin(c1*dt*wl)./wl.*v0_hat;
    temp1(1) = c1*dt*v0_hat(1) + u0_hat(1);
    u1_hat = temp1;
    u1 = ifft(temp1);
    
    temp2 = -u0_hat.*wl.*sin(c1*dt*wl)+cos(c1*dt*wl).*v0_hat;
    temp2(1) = v0_hat(1);
    v1 = ifft(temp2);

    % ---------------------------------------------------------------------
    % step two, full timestep, dt
    u2 = u1;
    DM2sp_u1 = ifft(abs(mu').^(2*sp).*u1_hat);
    DM2s_u1 = ifft(abs(mu').^(2*s).*u1_hat);
    v2 = v1 + c2*dt*(alpha*(-DM2sp_u1+DM2s_u1)+lambda*u1.^3);
    
    % ---------------------------------------------------------------------
    % step three, half timestep, 0.5*dt
    u2_hat = fft(u2);
    v2_hat = fft(v2);
    
    temp1 = u2_hat.*cos(c1*dt*wl)+sin(c1*dt*wl)./wl.*v2_hat;
    temp1(1) = c1*dt*v2_hat(1) + u2_hat(1);
    u3_hat = temp1;
    
    temp2 = -u2_hat.*wl.*sin(c1*dt*wl)+cos(c1*dt*wl).*v2_hat;
    temp2(1) = v2_hat(1);
    v3_hat = temp2;
    
    
    u0_hat = u3_hat;
    v0_hat = v3_hat;
    
    tn = n*dt;
end
u = real(ifft(u3_hat));
save("data_s1_tau"+num2str(NT)+"_t"+num2str(tn)+".mat",'x','tn','u');