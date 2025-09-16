function [u, i, rs] = CG_Solver(alpha, dt, u_n, f, M, mu, s, sp, lambda)
% use CG method to solve linear system created by Crank-Nicolson scheme
% 
% u^{n+1} - 0.5*alpha*dt^2*D*u^{n+1} -0.5*f(u^{n+1}) 
%        = 0.5*dt^2*(alpha D*u^{n-1}+f(u^{n-1})) - u^{n-1} + 2u^{n}
% 
% with D*u computed with the accelerated approach
% 2022/08/10


% initial guess
u0 = u_n;

u0_hat = fft(u0);
if length(s) == 1
    LMu0 = ifft(abs(mu).^(2*s).*u0_hat);
else
    LMu0 = LMu(u0_hat, M, mu, s, sp);
end
% residual
r = f - (u0 - 0.5*dt^2*(alpha*LMu0+lambda*u0.^3));
p = r;
rsold = r'*r;

u = u0;
tol = 1;
i = 0;
while tol >= 5e-13
    p_hat = fft(p);
    if length(s) == 1
        LMp = ifft(abs(mu).^(2*s).*p_hat);
    else
        LMp = LMu(p_hat, M, mu, s, sp);
    end
    Ap = p - 0.5*dt^2*(alpha*LMp+lambda*u.^2.*p);
   
    beta = rsold/(p'*Ap);
    u = u + beta*p;
    % tol = sqrt(rsnew);
    % real residual as tol
    u_hat = fft(u);
    if length(s) == 1
        LMuu = ifft(abs(mu).^(2*s).*u_hat);
    else
        LMuu = LMu(u_hat, M, mu, s, sp);
    end
    r = f - (u - 0.5*dt^2*(alpha*LMuu+lambda*u.^3));
    rsnew = r'*r;
    tol = sqrt(rsnew);
    
    i = i + 1;
    p = r + (rsnew/rsold)*p;
    rsold = rsnew;
end

% -------------------------------------------------------------------------

u_hat = fft(u);
LMuu = LMu(u_hat, M, mu, s, sp);
% real residual
r = f - (u - 0.5*dt^2*(alpha*LMuu+lambda*u.^3));
rs = sqrt(r'*r);