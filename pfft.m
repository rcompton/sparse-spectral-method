function zhatdown = pfft(z,picks,bigN)
%same action as fft(R'*z)/sqrt(n)

[m n]= size(z);
n = max(m,n);

zhatdown = zeros(bigN,1);
zhatdown(picks) = z;

zhatdown = fft(zhatdown);
zhatdown = (1/sqrt(n))*zhatdown;

end