function error_ref(Nref)

a = 32;
AX = -a;
BX = a;
dts = [1/128 1/256 1/512 1/1024 1/2048 1e-4];
%Nref = 4096;
h = (BX-AX)/Nref;

ref = load("data_sin_tau10000_t1.mat");
fid = fopen("TS2_Temporal_N"+num2str(Nref)+"_sin.txt",'w');
for n = 1:length(dts)
    dt = dts(n);
    sol = load("data_sin_tau"+num2str(1/dt)+"_t1.mat");
    e_Inf = max(abs(ref.u - sol.u));
    e_L2 = sqrt(h*sum((ref.u - sol.u).^2));
    fprintf(fid, '%4.16f %4.32f %4.32f\n', dt, e_Inf, e_L2);
end
calculate_convergence_rate("TS2_Temporal_N"+num2str(Nref)+"_sin.txt")
        