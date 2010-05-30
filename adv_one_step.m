function out = adv_one_step(r,psi, V, dt, M, Lnot)


N = size(r);

N = N(1);


k = N*linspace(-1/2,1/2,N); k = k';

GK = fftshift(exp(-(1i*dt/(4*M))*((2*pi/Lnot)^2)*(k.^2))); %% dt/2 kinetic energy propagator
GK2 = fftshift(exp(-(1i*dt/(2*M))*((2*pi/Lnot)^2)*(k.^2))); %% dt kinetic energy propagator

GV = exp(-1i*dt*V); %% Potential spatial interaction

%
% Apply exp(-c*laplacian)
%
ipsi = fft(psi);
psi = ifft(ipsi.*GK2);

%
% Apply exp(-c*V)
%
psi = GV.*psi;

%
% Apply exp(-c*laplacian)
%
ipsi = fft(psi);
out = ifft(ipsi.*GK2);

