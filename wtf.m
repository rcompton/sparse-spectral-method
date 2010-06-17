%what the hell I swear this was working!
clear all;close all;clc;

N = 4258;
uHat = randn(1,N) + 1i*randn(1,N);
zs = randsample(1:N,.90*N);
uHat(zs) = zeros(size(zs));

uHatBird = ifft(uHat);

uHat_for_real = fft(uHatBird);

fprintf('Fourier transform error %f \n', norm(uHat - uHat_for_real));

%time to party
practice_points = sort(randsample(1:N, N/2.5));
uder = zeros(1,N);
uder(practice_points) = uHatBird(practice_points);

% FPC_AS A_operator class
% weird shit here...
% FPC_AS likes A*A^t but fft is unitary
% so we'd like to work with A*conj(A) instead...
% shit
A = A_operator( @(z) pifft(z, find(uder)), @(z) pfft(z, find(uder), N) );
mu = 1e-10;
[ue, ~] = FPC_AS(N, A, nonzeros(uder), mu);
ue = ue'*sqrt(N);
ue = conj(ue); %???!!!!

figure(1)
plot(real(uHat));hold on;plot(real(ue),'r');
figure(2)
plot(imag(uHat));hold on;plot(imag(ue),'r');

fprintf('FPC error: %f \n', norm(uHat - ue));

% %norm(uHat' - uHatder*sqrt(N))
% 
% %Hack to get the scale right,
% %Fourier transform is an isometry of L2
% %Pe = Pe/norm(Pe,2);
% %Pe = Pe*norm(Pt);
% 
% 
% % Hurr durr scale by T
% Pe = Pe*sqrt(T*pi);
% 
% %This...
% Pe = fftshift(Pe);
% 
% 
% 
% 
% samps = sort(randsample(1:N,.35*N));
% uHatBird_samps = zeros(1,N);
% uHatBird_samps(samps) = uHatBird(samps);
% % FPC_AS A_operator class
% A = A_operator( @(z) pifft(z,samps), @(z) pfft(z,samps,N) );
% mu = 1e-10;
% [uHatder, ~] = FPC_AS(N, A, nonzeros(uHatBird_samps), mu);
% 
% 
% norm(uHat' - uHatder*sqrt(N))