function Pe = fft_via_fpca(Po, sample_idxs)
%
% Do an fft using FPCA
% The input data is large and sparse
% (ie include the 0s)
%
n = max(size(Po));
if ~all(size(Po) == [n 1])
    Po = Po';
end
assert(all(size(Po) == [n 1]));

% FPC_AS A_operator class
A = A_operator( @(z) pifft(z,sample_idxs), @(z) pfft(z,sample_idxs,n) );
u_samples = Po(sample_idxs);
mu = 1e-10;
[Pe, ~] = FPC_AS(n, A, u_samples, mu);
