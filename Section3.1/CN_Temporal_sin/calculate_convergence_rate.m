function calculate_convergence_rate(filename)
% filename : filename.dat

% import the error data
A = importdata(filename);
[m,n] = size(A);

% compute the convergence rate
AA = zeros(m,n+2);
AA(:,1:n) = A(:,1:n);
AA(:,n+1) = A(:,n);
AA(2:end,n) = log(A(1:end-1,n-1)./A(2:end,n-1))./log(A(1:end-1,n-2)./A(2:end,n-2));
AA(2:end,n+2) = log(A(1:end-1,n)./A(2:end,n))./log(A(1:end-1,n-2)./A(2:end,n-2));

% save the error and convergence rate information
save(filename,'AA','-ASCII');
