function splitting2d(xb,xe,yb,ye,xx,yy,kappa,T0,T1,N,dt,s,sp,phi,psi,lambda,M)
% -------------------------------------------------------------------------
% SPLITTING2D Solving the 2D variable-order fractional wave equation with
%            timesplitting for temporal discretization
%
%   Input:
%       xb, xe   - Left and right spatial boundaries
%       yb, ye   - Lower and upper spatial boundaries
%       kappa    - Diffusion coefficient
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
%       September 15 2025 (Monday)
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
w_lm = sqrt(kappa)*mu_lm.^(sp);
datetime(now,'ConvertFrom','datenum')


u0 = phi;
v0 = psi;
save("data_coscos_N"+num2str(N)+"_t"+num2str(T0) ...
     +".mat",'T0','xx','yy','u0','N');
 
u0_hat = fft2(u0);
v0_hat = fft2(v0);
% matrix free version -----------------------------------------------------
s0 = 0.5*(max(max(s))+min(min(s)));

c1 = 0.5;
c2 = 1;
%figure
%view(2)
for n = 1:NT
    % ---------------------------------------------------------------------
    % step one, half timestep, 0.5*dt
    temp1 = u0_hat.*cos(c1*dt*w_lm)+sin(c1*dt*w_lm)./w_lm.*v0_hat;
    temp1(1) = c1*dt*v0_hat(1) + u0_hat(1);
    u1_hat = temp1;
    u1 = ifft2(temp1);
    
    temp2 = -u0_hat.*w_lm.*sin(c1*dt*w_lm)+cos(c1*dt*w_lm).*v0_hat;
    temp2(1) = v0_hat(1);
    v1 = ifft2(temp2);

    % ---------------------------------------------------------------------
    % step two, full timestep, dt
    u2 = u1;
    DM2sp_u1 = ifft2(abs(mu_lm).^(2*sp).*u1_hat);
    DM2s_u1 = LMu(u1_hat, M, mu_lm, s, s0);
    v2 = v1 + c2*dt*(-kappa*(-DM2sp_u1+DM2s_u1)+lambda*u1.^3);
    
    % ---------------------------------------------------------------------
    % step three, half timestep, 0.5*dt
    u2_hat = fft2(u2);
    v2_hat = fft2(v2);
    
    temp1 = u2_hat.*cos(c1*dt*w_lm)+sin(c1*dt*w_lm)./w_lm.*v2_hat;
    temp1(1) = c1*dt*v2_hat(1) + u2_hat(1);
    u3_hat = temp1;
    
    temp2 = -u2_hat.*w_lm.*sin(c1*dt*w_lm)+cos(c1*dt*w_lm).*v2_hat;
    temp2(1) = v2_hat(1);
    v3_hat = temp2;
    
    
    u0_hat = u3_hat;
    v0_hat = v3_hat;
    
    tn = n*dt;
end
u = real(ifft2(u3_hat));
save("data_coscos_N"+num2str(N)+"_t"+num2str(tn)+".mat", ...
             'xx','yy','tn','u');
end