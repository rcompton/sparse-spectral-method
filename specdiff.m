function [ DDf ] = specdiff( f, x )
%spectral laplacian a column vector

n = length(x);
N = 2*n;
f = [f ;zeros(size(f))];

eLL = 2*(max(x) - min(x));
k = -N/2:(N/2-1);
k = k*2*pi/(eLL);
k = fftshift(k)';


DDf = ifft(-(k.^2).*fft(f));
DDf = DDf(1:n);

end

