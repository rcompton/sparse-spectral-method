%
% Fourier method solution to finite square well
%
clear all; close all;clc;

M = 1/2; %mass in au?

N = 512;

x = linspace(0.1, 50, N);
x = x';
dx = x(43)-x(42);

re = 3.28892;
De = 0.01688;
betah = 1.47612;

k = N*linspace(-1/2,1/2,N); k = k';
k = -N/2:N/2-1;
k = k';
k = k*1*pi/abs(max(x)-min(x));
k = fftshift(k);

%Working dimensions
dt = 1e-4 %% Time step
NPt = 2^15
T = dt*NPt

V = De*(1 - exp(-betah*(x-3*re))).^2 - De;
V = V - V.*(V>1.5) + 1.5.*(V>0.5);

V = -(1:length(x) < 200 ).*( 1:length(x) > 75);
V = 50*V';


Phi0 = exp(-(((x-2*re)/0.25).^2)); 
Phi0 = Phi0/trapz(x,Phi0);


%Uh...
%k = (pi/a)*(-N/2:N/2-1); % pi/L is dk
%k = fftshift(k)';

%use Strang splitting to move to Psi1
GKfast = exp(-(1i*dt/(4*M))*(1)*(k.^2)); %% dt/2 kinetic energy propagator
GVfast = exp(-1i*dt*V); %% Potential spatial interaction

% momentum space propagation
iPhi1 = fft(Phi0).*GKfast;
% move into physical space and apply potential operator
Phi1 = ifft(iPhi1);
Phi1 = GVfast.*Phi1;
% move into momentum space and propagate kinetic energy
iPhi1 = fft(Phi1);
iPhi1 = iPhi1.*GKfast;
% move back into physical space
Phi1 = ifft(iPhi1);

Phis = zeros(length(Phi0), NPt);
Phis(:,1) = Phi0;
Phis(:,2) = Phi1;

%% Use the SOD to propagate the Fourier method
%H = inline( 'ifft(-(k.^2).*fft(Phi)) + V.*Phi','V','k','Phi'); 
H = inline( '1*del2(Phi,dx) + V.*Phi', 'V', 'dx', 'Phi');

tic
for n=3:NPt-1
    %Compute HPhinm1 via Fourier method
    %HPhinm1 =ifft(-(k.^2).*fft(Phis(:,n-1))) + V.*Phis(:,n-1);
    %HPhinm1 = H(V,k,Phis(:,n-1));
    HPhinm1 = H(V,dx,Phis(:,n-1));

    Phis(:,n) = Phis(:,n-2) - 2i * dt * HPhinm1;
    
%     % momentum space propagation
%     iPhi = fft(Phis(:,n)).*GKfast;
%     % move into physical space and apply potential operator
%     Phis(:,n) = ifft(iPhi);
%     Phis(:,n) = GVfast.*Phis(:,n);
%     % move into momentum space and propagate kinetic energy
%     iPhi = fft(Phis(:,n));
%     iPhi = iPhi.*GKfast;
%     % move back into physical space
%     Phis(:,n+1) = ifft(iPhi);
    
    
end
toc

imshow(real( Phis(:,1:3:max(size(Phis))/3 ))) ;
colormap hot;

%% Make the autocorrelation function
Pt = zeros(1,NPt);
for nrn=1:NPt
    Pt(nrn) = trapz(x,Phis(:,nrn).*conj(Phi0));
end

t = linspace(0,T,NPt);
E = (1/dt)*(linspace(-pi,pi,length(Pt)));
Pt = (1-cos(2*pi*t/T)).*Pt;

Pe = fft(Pt);
Pe = fftshift(Pe)/T;

figure();
plot(E,abs(Pe));

figure();
plot(x,V,'rx');hold on;
plot(x,real(Phis(:, 1:72:1550)));
