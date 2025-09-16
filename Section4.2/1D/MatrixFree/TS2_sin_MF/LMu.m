function result = LMu(u_hat, M, mu, s, sp)

% assume all are column vectors
result = 0;

% -------------------------------------------------------------------------
% contirbution w.r.t. l=0 is zero, i.e., the first row of coe1 matrix is
% zero
coe(:,1) = abs(mu).^(2*sp);
for m = 1:M
    mm = [1:m];
    vec = (log(mu.^2)).^ones(1,m)./mm;
    coe(:,m+1) = prod(vec,2).*abs(mu).^(2*sp); % producation of each row
end
coe(1,:) = zeros(1,M+1);


for m = 0:M
    result = result + (s-sp).^m.*ifft(coe(:,m+1).*u_hat);
end