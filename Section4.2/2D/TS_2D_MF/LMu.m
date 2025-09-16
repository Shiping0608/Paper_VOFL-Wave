function result = LMu(u_hat, Q, mu_lm, s, sp)

[K,J] = size(mu_lm);
for k = 1:K
    for j = 1:J
        if mu_lm(k,j) == 0
            mu_lm(k,j) = 1;
        end
    end
end

coe = zeros(K,J,Q+1);
coe(:,:,1) = abs(mu_lm).^(2*sp);
for q = 1:Q
    for k = 1:q
        qq(1,1,k) = k;
    end
    vec = (log(mu_lm.^2)).^ones(1,1,q)./qq;
    coe(:,:,q+1) = prod(vec,3).*abs(mu_lm).^(2*sp);
end
%%%%

result = 0;
for q = 0:Q
    result = result + (s-sp).^q .* ifft2(coe(:,:,q+1).*u_hat);
end