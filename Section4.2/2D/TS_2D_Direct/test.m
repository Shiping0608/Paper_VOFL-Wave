function test

a = 32;
xb = -a;
xe = a;
yb = -a;
ye = a;
kappa = 0.2;
T0 = 0;
T1 = 1;

dt = 1e-4;
NN = [4 128 256 512 1024 2048];
for n = 1:length(NN)
    N = NN(n);
    h = (xe-xb)/N;
    x = [xb:h:xe-h];
    y = [yb:h:ye-h];
    [xx,yy] = meshgrid(x,y);
    % ---------------------------------------------------------------------
    s = 1-0.4*cos(0.25*pi*xx).*cos(0.25*pi*yy);
    sp = 0.5*(max(max(s))+min(min(s)));
    
    % ---------------------------------------------------------------------
    % initial condition phi and psi
    phi = 3*exp(-5*(xx.^2+yy.^2));
    psi = zeros(size(xx));
    
    %
    lambda = 0; % linear case
    diary("TS2_time_N"+num2str(N)+".txt");
    fid = fopen("TS2_time_N"+num2str(N)+".txt",'w');
    [t1, t2, t3] = fractionalWave_splitting2d(xb,xe,yb,ye,xx,yy,kappa, ...
                        T0,T1,N,dt,s,sp,phi,psi,lambda);
    fprintf(fid, '%6.10f %6.10f %6.10f\n', t1, t2, t3);
    fclose(fid);
end   
% error_ref()
