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
% all inputs are column vector, except M and sp

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