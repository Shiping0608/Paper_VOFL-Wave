function test()

a = 32;
xb = -a;
xe = a;
alpha = -1;
T0 = 0;
T1 = 1;

dt = 1e-4;
NN = [128 256 512 1024 2048 4096 8192];

for n = 1:length(NN)
    N = NN(n)
    
    diary("TS2_time_N"+num2str(N)+".txt");
    fid = fopen("TS2_time_N"+num2str(N)+".txt",'w');
    [t1, t2, t3] = fractionalWave_splitting(xb, xe, N, T0, T1, dt, alpha);
    fprintf(fid, '%6.10f %6.10f %6.10f\n', t1, t2, t3);
    fclose(fid);
end
