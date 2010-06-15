function Pe = fft_via_fpca(Po)
%
% Do a sparse fft using FPC
% The input data is large and sparse
% (ie include the 0s)
%
n = max(size(Po));
tflag = 0;
if ~all(size(Po) == [n 1])
    Po = Po';
    tflag = 1;
end
assert(all(size(Po) == [n 1]));
assert(nnz(Po) < .5*n); % otherwise use an fft

sample_idxs = find(Po);

% FPC_AS A_operator class
A = A_operator( @(z) pifft(z,sample_idxs), @(z) pfft(z,sample_idxs,n) );
u_samples = Po(sample_idxs);
mu = 1e-10;
[Pe, ~] = FPC_AS(n, A, u_samples, mu);

if tflag
    Pe = Pe';
end