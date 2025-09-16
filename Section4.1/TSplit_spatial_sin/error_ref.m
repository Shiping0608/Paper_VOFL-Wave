function error_ref()

a = 32;
AX = -a;
BX = a;
NN = [64 128 256 512];

Nref = 4096;
ref = load("data_sin_N"+num2str(Nref)+"_t1.mat");
fid = fopen("TS2_Spatial_N"+num2str(Nref)+"_sin.txt",'w');
for n = 1:length(NN)
    N = NN(n);
    h = (BX-AX)/N;
    sol = load("data_sin_N"+num2str(N)+"_t1.mat");
    e_Inf = max(abs(ref.u(1:2^(7-n):end) - sol.u));
    e_L2 = sqrt(h*sum((ref.u(1:2^(7-n):end) - sol.u).^2));
    fprintf(fid, '%4.16f %4.32f %4.32f\n', h, e_Inf, e_L2);
end
calculate_convergence_rate("TS2_Spatial_N"+num2str(Nref)+"_sin.txt")
        