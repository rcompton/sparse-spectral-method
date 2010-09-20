close all;clear all;clc;
% Apply e^At

n = 170;
A = randn(n) + 1i*randn(n);
H = A*A';


% we want to apply e^(-iHdt) to psi0
dt = .10;
psi0 = 1.2*abs(randn(n,1)).*sin(linspace(-pi,pi,n))';
psi0 = psi0/norm(psi0);

%Chebyshev expansion

%normalize H to have unit spectral radius
%Wtf, original paper had a typo for this?!
%is this just a cultural thing?
eigh = eig(H);
R = dt*(max(eigh)-min(eigh))/2; %What's the 2 for?
G = min(eigh)*dt;

dE = max(eigh) - min(eigh);

%normalize H
%Hnorm = 2*(H - 0.5*(emax+emin)*eye(n))/(emax-emin);
%Hnorm = 2*(H - eye(n)*min(eigh))/(max(eigh) - min(eigh)) - eye(n);
%Hnorm = -1i*H*dt/R;

Hnorm = (2/dE)*H - ((2*min(eigh)/dE) + 1)*eye(n);

eig(Hnorm)

%%
%figure maxk for a given timestep
maxk = 20;
nexttenjays = 1:.3:1.5;
while max( abs(besselj(maxk,dE*dt*nexttenjays)) ) > 1e-6 %&& abs(besselj(maxk,dt*(jdt+pi))) >1e-6
    maxk = maxk+1;
    
end

maxk

%make maxk Tkpsis
Tkpsis = zeros(n,maxk);
Tkpsis(:,1) = psi0;
Tkpsis(:,2) = Hnorm*psi0;
for k=3:maxk
    Tkpsis(:,k) = 2*Hnorm*Tkpsis(:,k-1) - Tkpsis(:,k-2);
end

%make the sum
chebsum = besselj(0,dt*dE)*Tkpsis(:,1);
chebsum = chebsum + -2i*besselj(1,dt*dE)*Tkpsis(:,2)
for k=2:maxk-1
    chebsum = chebsum + 2*((-1i)^k)*besselj(k,dt*dE)*Tkpsis(:,k+1);
    bsls(k) = abs(besselj(k,dt*dE));
end
plot(bsls)


figure();
plot(real(expm(-1i*Hnorm*dE*dt)*psi0));
hold on;
plot(real(chebsum),'r--');