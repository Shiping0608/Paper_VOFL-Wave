function test_energy()

a = 32;
xb = -a;
xe = a;
alpha = -1;
T0 = 0;
T1 = 5;

dt = 1e-4;
NN = [4096];

M = 30;
for N = NN
    h = (xe-xb)/N;
    x = [xb:h:xe-h]';

    % ---------------------------------------------------------------------
    % s = 0.3
    s = 1+0.3*sin(0.125*pi*x);
    sp = 0.5*(max(s)+min(s));

    % ---------------------------------------------------------------------
    % initial conditions
    phi = exp(-x.^2);
    psi = zeros(size(x));
    
    splitting(xb, xe, N, T0, T1, dt, alpha, phi, psi, sp, s, M)
end


data = load("data_sin_N"+num2str(NN)+"_t3.mat");
t = [0:0.01:3]';
en_fd = data.Energy_FD(1:length(t));


% energy comparison
data1 = load("./../TSplit_spatial_s1/data_s1_N"+num2str(NN)+"_t3.mat");
en_fd1 = data1.Energy_FD(1:length(t));
figure
plot(t,en_fd1,'r-','linewidth',1.5,'DisplayName','$\bar{s}_{1}=1$'); hold on
plot(t,en_fd,'b-.','linewidth',1.5,'DisplayName','$s_{1}(x)$')
legend('fontsize',14,'interpreter','latex')
xlabel('t','fontsize',16)
ylabel('E(t)','fontsize',16)
axis([0 3 0.4 0.45])