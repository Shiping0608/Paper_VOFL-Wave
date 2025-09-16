function result = LMu(u_hat, M, mu, s, sp)
% -------------------------------------------------------------------------
% LMu() Using the fast algorithm to evalute the variable-order fractional
%       Laplacian
%
%   Input:
%       u_hat    - fft(uh)
%       M        - Number of modes / truncation parameter for LMu
%       mu       - DFT coefficients of the operator
%       T0, T1   - Initial and final time
%       s, sp    - Vector for s(x) and its average value
%
%   Output:
%       variable-order fractional Laplacian
%
%   Author:
%       Shiping Zhou
%
%   Date last updated:
%       [2025-09-15]
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
coe = (log(mu.^2)).^mm./factorial(mm).*abs(mu).^(2*sp);
coe(1,1) = 0;


for m = 0:M
    result = result + (s-sp).^m.*ifft(coe(:,m+1).*u_hat);
end