%% You can do something like an FFT using L1
%
clear all;
close all;
clc;

n = 2^11;
sparsity = 0.075*n;
num_samples = round(n/3.7321);

%seed the random generator with the example seeder
stream = RandStream('mrg32k3a');

%sparse freqs
uHat_exact = zeros(n,1);
target_points = randsample(stream,1:n,sparsity);
uHat_exact(target_points) = randn(stream,sparsity,1)*10;

%full time data (never used)
u = ifft(uHat_exact);

%make our downsampling matrix

sample_points = randsample(stream,1:n, num_samples);
R = zeros(num_samples, n);
for i=1:num_samples
    R(i,sample_points(i)) = 1;
end

%% reconstruct from sparse data

% FPC_AS A_operator class
A = A_operator( @(z) (1/sqrt(n))*R*fft(z), @(z) sqrt(n)*ifft(R'*z));

% tiny mu corresponds to heavy weight on the fidelity term
mu = 1e-10;
u_samples = A*uHat_exact; %apply the first operator in A_operator.

%u_tester = n*R*ifft(uHat_exact);
%norm(u_tester - u_samples)

% Call Wotao's code.
[uHat_approx, Out] = FPC_AS(n,A,u_samples,mu);


%%
figure(1)
title('exact Blue, approx red')
plot(1:n, uHat_exact,'b')
hold on
plot(1:n, uHat_approx,'r')

fprintf('L2 error in approx: %f \n', norm(uHat_approx - uHat_exact,2));%/norm(uHat_exact));
