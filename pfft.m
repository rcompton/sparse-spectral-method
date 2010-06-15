function zhatdown = pfft(z,picks,bigN)
%same action as fft(R'*z)/sqrt(n)

n = max(size(z));

zhatdown = zeros(bigN,1);

%size(z)
%size(zhatdown(picks))

zhatdown(picks) = z;

zhatdown = fft(zhatdown);
zhatdown = (1/sqrt(n))*zhatdown;

end