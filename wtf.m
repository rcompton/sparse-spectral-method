
%what the hell I swear this was working!

N = 24258;
uHat = randn(1,N);
zs = randsample(1:N,.90*N);
uHat(zs) = zeros(size(zs));

uHatBird = ifft(uHat);

uHat_for_real = fft(uHatBird);

norm(uHat - uHat_for_real)

%time to party
samps = sort(randsample(1:N,.35*N));
uHatBird_samps = zeros(1,N);
uHatBird_samps(samps) = uHatBird(samps);
% FPC_AS A_operator class
A = A_operator( @(z) pifft(z,samps), @(z) pfft(z,samps,N) );
mu = 1e-10;
[uHatder, ~] = FPC_AS(N, A, nonzeros(uHatBird_samps), mu);


norm(uHat' - uHatder*sqrt(N))
