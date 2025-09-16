function [t1, t2, t3] = splitting(xb, xe, N, T0, T1, dt, alpha, M)

tic

c1 = 0.5;
c2 = 1;
lambda = 1;

NT = round((T1-T0)/dt); % number of time levels
h = (xe-xb)/N;
x = [xb:h:xe-h]; % mesh

% ---------------------------------------------------------------------
s = 1+0.3*sin(0.125*pi*x)';
sp = 0.5*(max(s)+min(s));
s0 = 0.5*(max(s)+min(s));

% ---------------------------------------------------------------------
% initial conditions
u0 = exp(-x.^2)';
v0 = zeros(size(x))';
u0_hat = fft(u0);
v0_hat = fft(v0);


mu = [0:N/2-1,-N/2:-1]*(2*pi/(xe-xb));
wl = sqrt(-alpha)*abs(mu').^(sp(1));

t1 = toc;

% timing section 2 time for solution at dt --------------------------------
tic
t2 = toc;

% timing section 3 time for solutions at rest time levels -----------------
tic
for n=1:NT
    % ---------------------------------------------------------------------
    % step one, half timestep, 0.5*dt
    temp1 = u0_hat.*cos(c1*dt*wl)+sin(c1*dt*wl)./wl.*v0_hat;
    temp1(1) = c1*dt*v0_hat(1) + u0_hat(1);
    u1 = ifft(temp1);
    
    temp2 = -u0_hat.*wl.*sin(c1*dt*wl)+cos(c1*dt*wl).*v0_hat;
    temp2(1) = v0_hat(1);
    v1 = ifft(temp2);

    % ---------------------------------------------------------------------
    % step two, full timestep, dt
    %u2 = u1;
    % v2 = v1 + dt*alpha*(-DM2_sp+DM2_s)*u1;
    % accelerated version 8/7/2022
    DM2sp_u1 = ifft(abs(mu').^(2*sp).*temp1);
    DM2s_u1 = LMu(temp1, M, double(mu'), double(s), s0);
    v2 = v1 + c2*dt*(alpha*(-DM2sp_u1+DM2s_u1)+lambda*u1.^3);
    
    % ---------------------------------------------------------------------
    % step three, half timestep, 0.5*dt
    u2_hat = fft(u1);
    v2_hat = fft(v2);
    
    temp1 = u2_hat.*cos(c1*dt*wl)+sin(c1*dt*wl)./wl.*v2_hat;
    temp1(1) = c1*dt*v2_hat(1) + u2_hat(1);
    
    temp2 = -u2_hat.*wl.*sin(c1*dt*wl)+cos(c1*dt*wl).*v2_hat;
    temp2(1) = v2_hat(1);
    
    u0_hat = temp1;
    v0_hat = temp2;
end
t3 = toc;