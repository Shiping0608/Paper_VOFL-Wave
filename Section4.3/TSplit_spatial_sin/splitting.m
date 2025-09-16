function splitting(xb, xe, N, T0, T1, dt, alpha, phi, psi, sp, s, M)


NT = round((T1-T0)/dt); % number of time levels
h = (xe-xb)/N;
x = [xb:h:xe-h]; % mesh

mu = [0:N/2,-N/2+1:-1]*(2*pi/(xe-xb));
wl = sqrt(-alpha)*abs(mu').^(sp(1));

u0 = phi; % u(x,0)
v0 = psi; % u_{t}(x,0)

u0_hat = fft(u0);
v0_hat = fft(v0);

s0 = 0.5*(max(s)+min(s));

c1 = 0.5;
c2 = 1;
lambda = 1;

dtt = 0.01;
Ndtt = round((T1-T0)/dtt);

e0 = @(x) 2*x.^2.*exp(-2*x.^2)-1/4*exp(-4*x.^2);
Energy_FD = zeros(Ndtt+1,1);
Energy_SP = zeros(Ndtt+1,1);
Energy_FD(1) = integral(e0,xb,xe); % same value with Trapezoidal rule
Energy_SP(1) = integral(e0,xb,xe);

unm1 = u0;

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
    DM2s_u1 = LMu(u1_hat, M, double(mu'), double(s), s0);
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
    if tn == 0.01
        un = real(ifft(u3_hat));
    end
    
    if tn >= 0.02 && mod(tn,0.01) == 0
        nk = round(tn/0.01);
        unp1 = real(ifft(u3_hat));
        u_t = abs((unp1-unm1)/2/dtt);
        u_x1 = abs(([un(2:N);un(1)]-[un(N);un(1:N-1)])/2/h);
        u_x2 = abs(ifft(1i*mu'.*fft(un)));
        Energy_FD(round(nk)) = ...
            h/2*(1/2*u_t(1)^2 + 1/2*u_x1(1)^2 - 1/4*un(1)^2) ...
          + h*sum(1/2*u_t(2:N-1).^2 + 1/2*u_x1(2:N-1).^2 - 1/4*un(2:N-1).^4) ...
          + h/2*(1/2*u_t(1)^2 + 1/2*u_x1(1)^2 - 1/4*un(1)^2);
        Energy_SP(round(nk)) = ...
            h/2*(1/2*u_t(1)^2 + 1/2*u_x2(1)^2 - 1/4*un(1)^2) ...
          + h*sum(1/2*u_t(2:N-1).^2 + 1/2*u_x2(2:N-1).^2 - 1/4*un(2:N-1).^4) ...
          + h/2*(1/2*u_t(1)^2 + 1/2*u_x2(1)^2 - 1/4*un(1)^2);

        unm1 = un;
        un = unp1;
    end
end
u = real(ifft(u3_hat));
save("data_sin_N"+num2str(N)+ ...
             "_t"+num2str(tn)+".mat",'x','tn','u','Energy_SP','Energy_FD');