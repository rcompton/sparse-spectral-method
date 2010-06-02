function zbirddown = pifft(z,picks)
%same action as R*ifft(z)*sqrt(n)

[m n]= size(z);
n = max(m,n);

zbirddown = ifft(z);
zbirddown = zbirddown(picks);
zbirddown = sqrt(n)*zbirddown;
end