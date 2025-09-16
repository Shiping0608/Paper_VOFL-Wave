function result = LMu(u_hat, M, mu, s, sp)

% assume all are column vectors
result = 0;

for k = 1:length(mu)
    if mu(k) == 0
        mu(k) = 1;
    end
end

% -------------------------------------------------------------------------
% first row of coe == 0, to ensure that mu_l has no contribution for l=0
% term
mm = [0:M];
coe = (log(mu.^2)).^mm./factorial(mm).*abs(mu).^(2+2*sp);
coe(1,1) = 0;


for m = 0:M
    result = result + (s-sp).^m.*ifft(coe(:,m+1).*u_hat);
end