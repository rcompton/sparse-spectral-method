% 
% It's called invFT, but I want an FT....
%
% In this example we know that the frequency space data
% is sparse and we get it from sparse physical space data
% using FPC_AS
%
clear all;
close all;
clc;

n = 2^9.14;
sparsity = 0.075*n;
num_samples = round(n/3.7321);

%seed the random generator with the example seeder
stream = RandStream('mrg32k3a');

%sparse freqs
uHat_exact = zeros(n,1);
target_points = randsample(stream,1:n,sparsity);
uHat_exact(target_points) = randn(stream,sparsity,1)*10;

%full time data (never used)
u = ifft(uHat_exact)*sqrt(n);

%make our downsampling matrix

sample_points = randsample(stream,1:n, num_samples);
R = zeros(num_samples, n);
for i=1:num_samples
    R(i,sample_points(i)) = 1;
end

%% reconstruct from sparse data
A = A_operator( @(z) pifft(z,sample_points), @(z) pfft(z,sample_points,n) );

% FPC_AS A_operator class
%A = A_operator( @(z) (sqrt(n))*R*ifft(z), @(z) (1/sqrt(n))*fft(R'*z));

% tiny mu corresponds to heavy weight on the fidelity term
mu = 1e-10;
%u_samples = A*uHat_exact; %apply the first operator in A_operator (an inverse fft)
u_samples = u(sample_points);

%u_tester = sqrt(n)*R*ifft(uHat_exact);
%norm(u_tester - u_samples)

% Call Wotao's code.
[uHat_approx, Out] = FPC_AS(n,A,u_samples,mu);


%%
figure(1)
title('exact Blue, approx red')
plot(1:n, uHat_exact,'b')
hold on
plot(1:n, uHat_approx,'r')

%The point is that uHat_approx is a fourier transform of u
fprintf('nonrelative L2 error in approx: %f \n', norm(uHat_approx - uHat_exact,2));%/norm(uHat_exact));
fprintf('L2 error a different way %f \n', norm(fft(u)/sqrt(n) - uHat_approx,2));

