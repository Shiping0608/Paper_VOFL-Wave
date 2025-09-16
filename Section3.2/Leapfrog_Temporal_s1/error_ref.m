function error_ref()

a = 32;
xb = -a;
xe = a;
dts = [1/128 1/256 1/512 1/1024 1/2048];

Nref = 4096;
h = (xe-xb)/Nref;
t= 1;
ref = load("sol_s0_tau10000_t"+num2str(t)+".mat");
fid = fopen("Temporal_Leapfrog_dt_ref0.0001" ...
            +"_s0_t"+num2str(t)+".txt",'w');
for n = 1:length(dts)
    dt = dts(n);
    sol = load("sol_s0_tau"+num2str(1/dt)+"_t"+num2str(t)+".mat");
    e_Inf = max(abs(ref.un - sol.un));
    e_L2 = sqrt(h*sum((ref.un - sol.un).^2));
    fprintf(fid, '%4.16f %4.32f %4.32f\n', dt, e_Inf, e_L2);
end
calculate_convergence_rate("Temporal_Leapfrog_dt_ref0.0001" ...
                           +"_s0_t"+num2str(t)+".txt")
        