%
% In this example we know that the frequency space data
% is sparse and we get it from sparse physical space data
% using FPC_AS
% without making an explicit sample matrix R
%
clear all;
close all;
clc;

n = 2^10.42;
sparsity = 0.0612*n;
num_samples = round(n/3.8121);

%seed the random generator with the example seeder
stream = RandStream('mrg32k3a');

%create the sparse freqs that we want to reconstruct
uHat_exact = zeros(n,1);
target_points = randsample(stream,1:n,sparsity);
uHat_exact(target_points) = randn(stream,sparsity,1)*10;


%full time data (never used)
u = sqrt(n)*ifft(uHat_exact);

%make our downsampling matrix
sample_points = randsample(stream,1:n, num_samples);

%% reconstruct from sparse data

% FPC_AS A_operator class
A = A_operator( @(z) pifft(z,sample_points), @(z) pfft(z,sample_points,n) );

% tiny mu corresponds to heavy weight on the fidelity term
mu = 1e-10;

%downsample.
%u_samples are what we measure from the expensive machine
%we simulate measured data with the following line:
u_samples = u(sample_points); %apply the first operator in A_operator.

% Call Wotao's code.
[uHat_approx, Out] = FPC_AS(n,A,u_samples,mu);


%%
figure(1)
title('exact Blue, approx red')
plot(1:n, uHat_exact,'b')
hold on
plot(1:n, uHat_approx,'r')

fprintf('L2 error in approx: %f \n', norm(uHat_approx - uHat_exact,2));%/norm(uHat_exact));
