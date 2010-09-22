close all;clear all;clc;
% Apply e^(lap + V) over a timestep dt

n = 512; % number of divisions of x (should be even!!)
L = 50;

x = linspace(.1, L, n);
dx = x(21) - x(20)

%steps = 200; % steps calculated per frame dt < .25*dx^2 ?Depends on M
%dt = 1 / steps; %time step, unstable at 1/200 with SOD!
dt = 73;

re = 3.2889;
De = 0.01688;
betah = 1.47612;

u = exp(-((x-3)/.25).^2);

V = De*(1-exp(-betah*(x-re))).^2 - De;
V = V.*(V < 0.5) + (V>0.5);
V = V';

% Wtf?
M = 50;

dE = (pi^2)/(2*M*dx^2) + max(V) - min(V)
maxE = (pi^2)/(2*M*dx^2) + max(V)
minE = min(V)

psi0 = exp(-((x-3)/.25).^2);
%psi0 = exp(-0.1*(x-10*L/11).^2);
psi0 = psi0';

%Chebyshev expansion

%normalize H to have unit spectral radius
%Wtf, original paper had a typo for this?!
%is this just a cultural thing?
%R = dt*(max(eigh)-min(eigh))/2; %What's the 2 for?
%G = min(eigh)*dt;

%Hnorm = (2/dE)*H - ((2*min(eigh)/dE) + 1)*eye(n);
%Hnorm = inline( '(2/dE)*((-1/(2*M))*4*del2(Phi,dx) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','dx','Phi');

eLL = max(x) - min(x);
k = -n/2:(n/2-1);
k = k*2*pi/(eLL);
k = fftshift(k)';

%Hnormspec = inline('(2/dE)*((-1/(2*M))* 4*del2(Phi,dx)          + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','k','dx','Phi');

Hnormspec = inline( '(2/dE)*((-1/(2*M))* specdiff(Phi,x) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','k','dx','x','Phi');
%Hnormspec = inline( '(2/dE)*((-1/(2*M))* ifft(-k.^2 .*fft(Phi)) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','k','dx','x','Phi');




%%

%figure maxk for a given timestep
maxk = 2;
nexttenjays = 1:.3:1.5;
while max( abs(besselj(maxk,dE*dt*nexttenjays)) ) > 1e-6 %&& abs(besselj(maxk,dt*(jdt+pi))) >1e-6
    maxk = maxk+1;    
end

maxk

%make maxk Tkpsis
Tkpsis = zeros(n,maxk);
Tkpsis(:,1) = psi0;
%Tkpsis(:,2) = Hnorm(M,V,dE,minE,dx,psi0);
Tkpsis(:,2) = Hnormspec(M,V,dE,minE,k,dx,x,psi0);
for kunt=3:maxk
    %Tkpsis(:,k) = 2*Hnorm(M,V,dE,minE,dx,Tkpsis(:,k-1)) - Tkpsis(:,k-2);
    Tkpsis(:,kunt) = 2*Hnormspec(M,V,dE,minE,k,dx,x,Tkpsis(:,kunt-1)) - Tkpsis(:,kunt-2);
end

%make the sum
chebsum = besselj(0,dt*dE)*Tkpsis(:,1);
chebsum = chebsum + -2i*besselj(1,dt*dE)*Tkpsis(:,2);
for kunt=2:maxk-1
    chebsum = chebsum + 2*((-1i)^kunt)*besselj(kunt,dt*dE)*Tkpsis(:,kunt+1);
    bsls(kunt) = abs(besselj(kunt,dt*dE));
end
figure()
plot(bsls)

%check it

figure();
%plot(real(expm(-1i*Hnorm*dE*dt)*psi0));
hold on;
plot(real(chebsum),'r--');
plot(imag(chebsum),'b--');
plot(abs(chebsum),'k','Linewidth',2);
plot(V,'g','Linewidth',3);
