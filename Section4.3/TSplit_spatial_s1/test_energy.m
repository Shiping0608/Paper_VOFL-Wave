function test_energy()

a = 32;
xb = -a;
xe = a;
alpha = -1;
T0 = 0;
T1 = 3;

dt = 1e-4;
NN = [4096];


for N = NN
    h = (xe-xb)/N;
    x = [xb:h:xe-h];

    % ---------------------------------------------------------------------
    s = 1;
    sp = 0.5*(max(s)+min(s));

    % ---------------------------------------------------------------------
    % initial conditions
    phi = exp(-x.^2)';
    psi = zeros(size(x))';
    
    splitting(xb, xe, N, T0, T1, dt, alpha, phi, psi, sp, s)
end
