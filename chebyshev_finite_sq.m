%
% Fourier-Chebyshev method solution to finite square well
%
clear all; close all;clc;
a = 48;  %% Size of domain
M = 1/2; %% Mass (always is 0.5?)
N = 512; %% Spatial domain

x = linspace(-a/2,a/2,N);
x = x';
dx = x(34)-x(33);

%k = N*linspace(-1/2,1/2,N); k = k';
k = -N/2:N/2-1;
k = k';
k = k*1*pi/a;
k = fftshift(k);

%Working dimensions
dt = 1e-3; %% Time step
NPt = 200000;
T = dt*NPt;

V0 = 200;
V = zeros(length(x),1) - V0; % 1*((2*x).^2 - (0.6*a)^2);  %  
b = a/16;
V(x<-b) = 0;
V(x>+b) = 0;
% V(1:5) = V(1:5)-i*1e3;  %% Absorption at simulation boundary
% V(end-5:end) = V(end-5:end)-i*1e3;  %% Absorption at simulation boundary

%% Get initial data...

Phi0 = exp(-(5*(x-0*a/128)).^2); 


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

%% Use the Chebyshev to propagate the Fourier method
%GK =
inline(sprintf('fftshift(exp(-(i*dt/(4*%f))*((2*pi/%f)^2)*(k.^2)))',M,a)','dt','k'); %% dt/2 kinetic energy propagator

H = inline( 'ifft(-(k.^2).*fft(Phi)) + V.*Phis','V','k','Phi'); 
    
tic
for n=3:NPt-1
    
    %Compute HPhinm1 via Fourier method
    HPhinm1 = ifft(-(k.^2).*fft(Phis(:,n-1))) + V.*Phis(:,n-1); 
    

    
    nterms = 50;
    
    
    
    %Phis(:,n) = Phis(:,n-2) - 2i * dt * HPhinm1;
    
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

imshow(real( Phis(:,1:5:max(size(Phis))/10 ))) ;
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

