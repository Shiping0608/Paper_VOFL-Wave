function result = LMu(u_hat, Q, mu_lm, s, sp)
% LMU   Apply fractional spectral operator in Fourier space.
%
%   result = LMu(u_hat, Q, mu_lm, s, sp)
%
%   This function evaluates a fractional operator in Fourier space,
%   expanded as a truncated logarithmic series.
%
%   Input:
%       u_hat  - Fourier coefficients of the function u
%       Q      - Maximum order of truncated logarithmic expansion
%       mu_lm  - Frequency matrix (K x J), e.g. from meshgrid of wavenumbers
%       s      - Variable-order
%       sp     - Averaged-order (also the choice for TS)
%
%   Output:
%       result - Variable-order fractional Laplacian
%
%
%   AUTHOR:
%       Shiping Zhou
%
%   DATE:
%       September 15, 2025 (Monday)
%
%   References:
%       - Y. Zhang, X. Zhao, and S. Zhou, "Fourier pseudospectral methods 
%         for the variable-order space fractional wave equations"
%         (under review, to be updated)
%
%   Licensing:
%    This code is distributed under the MIT license.
%
% -------------------------------------------------------------------------

[K,J] = size(mu_lm);
for k = 1:K
    for j = 1:J
        if mu_lm(k,j) == 0
            mu_lm(k,j) = 1;
        end
    end
end

%q = 0:Q;
coe1 = zeros(K,J,Q+1);
for q = 0:Q
    coe1(:,:,q+1) = (log(mu_lm.^2)).^q./factorial(q).*abs(mu_lm).^(2*sp);
end
%%%%

result = 0;
for q = 0:Q
    result = result + (s-sp).^q .* ifft2(coe1(:,:,q+1).*u_hat);
end