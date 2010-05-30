function out = psi_not(r)
%
% r goes from 0.1 to 50au
%
% propagate in time for Nt/4 points
% Nt = 2^13
% tau = 4

out = exp(-((r-3.0)./0.25).^2);

%wavefunctions have L2 mass 1.0
out = out.*sqrt(1.0/trapz(r,abs(out).^2));