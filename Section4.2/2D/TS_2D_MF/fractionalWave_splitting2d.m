function [t1, t2, t3] = fractionalWave_splitting2d(xb,xe,yb,ye,xx,yy, ...
                        kappa,T0,T1,N,dt,s,sp,phi,psi,lambda,M)

% solve the 2d fractional wave equation with central difference 
% time discretization with matrix-free acceleration
% 2022/06/1

s0 = 0.5*(max(max(s))+min(min(s)));

c1 = 0.5;
c2 = 1;

NT = round((T1-T0)/dt);

% timing section 1 --------------------------------------------------------
% differential matrix D2
tic
mu_l = [0:N/2-1,-N/2:-1]*(2*pi/(xe-xb));
mu_m = [0:N/2-1,-N/2:-1]*(2*pi/(ye-yb));
[mu_ll, mu_mm] = meshgrid(mu_l,mu_m);
mu_lm = sqrt(mu_ll.^2+mu_mm.^2);
w_lm = sqrt(kappa)*mu_lm.^(sp);
t1 = toc;

% timing section 2 --------------------------------------------------------
tic
u0 = phi;
v0 = psi;
 
u0_hat = fft2(u0);
v0_hat = fft2(v0);
t2 = toc;

% timing section 3 --------------------------------------------------------
tic
for n = 1:NT
    % ---------------------------------------------------------------------
    % step one, half timestep, 0.5*dt
    temp1 = u0_hat.*cos(c1*dt*w_lm)+sin(c1*dt*w_lm)./w_lm.*v0_hat;
    temp1(1) = c1*dt*v0_hat(1) + u0_hat(1);
    u1 = ifft2(temp1);
    
    temp2 = -u0_hat.*w_lm.*sin(c1*dt*w_lm)+cos(c1*dt*w_lm).*v0_hat;
    temp2(1) = v0_hat(1);
    v1 = ifft2(temp2);

    % ---------------------------------------------------------------------
    % step two, full timestep, dt
    u2 = u1;
    DM2sp_u1 = ifft2(abs(mu_lm).^(2*sp).*temp1);
    DM2s_u1 = LMu(temp1, M, mu_lm, s, s0);
    v2 = v1 + c2*dt*(-kappa*(-DM2sp_u1+DM2s_u1)+lambda*u1.^3);
    
    % ---------------------------------------------------------------------
    % step three, half timestep, 0.5*dt
    u2_hat = fft2(u2);
    v2_hat = fft2(v2);
    
    temp1 = u2_hat.*cos(c1*dt*w_lm)+sin(c1*dt*w_lm)./w_lm.*v2_hat;
    temp1(1) = c1*dt*v2_hat(1) + u2_hat(1);
    
    temp2 = -u2_hat.*w_lm.*sin(c1*dt*w_lm)+cos(c1*dt*w_lm).*v2_hat;
    temp2(1) = v2_hat(1);
    
    
    u0_hat = temp1;
    v0_hat = temp2;
end
t3 = toc;