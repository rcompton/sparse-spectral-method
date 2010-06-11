%
% In this example we know that the spectral data
% is sparse and we get it from sparse time data
% using FPC_AS
% without making an explicit sample matrix R
% and we use the Pt generated by copy_paste.m
%

n = Nt;
num_samples = round(n/17.8121);
num_sparsity = .05*n;

assert(any(size(Po)==Nt));

%seed the random generator with the example seeder
stream = RandStream('mrg32k3a');

%% Draw the energy spectrum with classical FFT methods
%Pe_old = fft(((1-cos(2*pi*t'/T)).*Po))/sqrt(n); %Hanning windows...
%Peww_old =  fft(Po)/sqrt(n); %no windowing in compressed sensing?


%holy crap....
%this was being reassigned on each run...
%things would have matrix mult errors every other time!
if ~all(size(Po) == [Nt 1])
    Po = Po';
end
assert(all(size(Po) == [Nt 1]));

%window Po
Po = (1-cos(2*pi*t'/T)).*Po;

%uHat_exact = zeros(n,1);
%target_points = randsample(stream,1:n,num_sparsity);
%uHat_exact(target_points) = randn(stream,num_sparsity,1)*10;

uHat_exact = fft(Po)/sqrt(n);

%full time data (never used)
%u = sqrt(n)*ifft(uHat_exact);

%% reconstruct from sparse data

%make our downsampling matrix
sample_points = randsample(stream,1:n, num_samples);

% FPC_AS A_operator class
A = A_operator( @(z) pifft(z,sample_points), @(z) pfft(z,sample_points,n) );

% tiny mu corresponds to heavy weight on the fidelity term
mu = 1e-10;

%downsample.
%u_samples are what we measure from the expensive machine
%we simulate measured data with the following line:
%u_samples = Po(sample_points)'; %apply the first operator in A_operator.
%u = ifft(uHat_exact)*sqrt(n);

%u_samples = u(sample_points);
%Po = Po';
u_samples = Po(sample_points);


%fprintf('u-Po %f\n',norm(u - Po));



% Call Wotao's code.
%[Pe_new, Out] = FPC_AS(n, A, u_samples, mu);

[uHat_approx, Out] = FPC_AS(n, A, u_samples, mu);

%Pe_new = Pe_new';

%%
%close all;
%figure(1);
%title('old Blue, new red');
%plot(1:n, abs(Peww_old*sqrt(NPt)),'b');
%hold on;
%plot(1:n, abs(Pe_new),'r--');
%hold on;
%plot(1:n, abs(Pe_old*sqrt(NPt)),'g');
%fprintf('L2 error in approx: %f \n', norm( abs(Pe_old*sqrt(NPt)) - abs(Pe_new),2)/norm(Pe_old*sqrt(NPt)) );

figure(1)
title('real parts')
plot(1:n, real(uHat_exact))
hold on
plot(1:n, real(uHat_approx),'r--')

figure(2)
title('imag parts')
plot(1:n, imag(uHat_exact))
hold on
plot(1:n, imag(uHat_approx),'r--')

fprintf('L2 unnormalized error %f \n', norm(uHat_approx - uHat_exact));
