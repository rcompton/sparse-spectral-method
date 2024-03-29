%% Finite Square Well: Split-Operator Method
%% Ref: M. D. Feit, J. A. Fleck, Jr., and A. Steiger, 
%%      "Solution of the Schrodinger Equation by a Spectral Method",
%%      Journal of Computational Physics 47, 412-433 (1982).
%% Ref: D. J. Griffiths, 
%%      "Introduction to Quantum Mechanics",
%%      ISBN 0-13-124405-1
%% 
clear all; close all;
a = 48;  %% Length
M = 1/2; %% Mass
N = 512;
x = linspace(-a/2,a/2,N); x = x';
%k = N*linspace(-1/2,1/2,N); k = k';
k = -N/2:N/2-1; k = k';
dt = 1e-3; %% Time step
%%-------------------------------------------------------------------------
%%
%%% Potential
V0 = 200;
V = zeros(length(x),1) - V0; % 1*((2*x).^2 - (0.6*a)^2);  %  
b = a/16;
V(x<-b) = 0;
V(x>+b) = 0;
% V(1:5) = V(1:5)-i*1e3;  %% Absorption at simulation boundary
% V(end-5:end) = V(end-5:end)-i*1e3;  %% Absorption at simulation boundary
%%%
%%-------------------------------------------------------------------------
%%
%%% initial wave packet
%% Gaussian
Phi0 = exp(-(5*(x-0*a/128)).^2);
%
%% Odd
% mode = 3;
% Phi0 = sqrt(2/b)*sin(mode*pi*x/b);
% Phi0(x<-b) = 0;
% Phi0(x>+b) = 0;
%
%% Even
% mode = 8;
% Phi0 = sqrt(2/b)*cos(mode*pi*x/b);
% Phi0(x<-b) = 0;
% Phi0(x>+b) = 0;
%%
%%-------------------------------------------------------------------------
%%
Phi0c = conj(Phi0); %% real(Phi0)- i*imag(Phi0);
%%
figure(1);set(gcf,'position',[37 208 538 732]);
plot(x,V,'r');hold on;plot(x,max(abs(real(V)))*abs(Phi0c));hold off; pause(1);
%%
GK = fftshift(exp(-(i*dt/(4*M))*((2*pi/a)^2)*(k.^2))); %% dt/2 kinetic energy propagator
GK2 = fftshift(exp(-(i*dt/(2*M))*((2*pi/a)^2)*(k.^2))); %% dt kinetic energy propagator
GV = exp(-i*dt*V); %% Potential spatial interaction
% plot((-(dt/(4*M))*((2*pi/a)^2)*(k.^2)));
% plot(-dt*V);
%%
iPhi = fft(Phi0);
Phi = ifft(iPhi.*GK);
Phi = GV.*Phi;
NPt = 50000;
Pt = zeros(1,NPt);
En = -105.99;  %% Energy eigen value
T = dt*NPt;
t = linspace(0,T,NPt);
wt = (1-cos(2*pi*t/length(t)));
uns = 0;
tic
for nrn = 1:NPt
    % momentum space propagation
    iPhi = iPhi.*GK;
    % move into physical space and apply potential operator
    Phi = ifft(iPhi);
    Phi = GV.*Phi;
    % move into momentum space and propagate again
    iPhi = fft(Phi);
    iPhi = iPhi.*GK;
    % move back into physical space and record P(t) at the sample point
    Phi = ifft(iPhi);
    Pt(nrn) = trapz(x, Phi0c.*Phi);
end
toc
iPhi = fft(Phi);
Phi = ifft(iPhi.*GK);
%%

%this is crap this is crap careful about what you're doing here!
% disp('!!!!!!!!warning shifting Pt!!');
% Pt = [Pt(1) Pt];
% Pt = Pt(1:length(Pt)-1);

estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
T = dt*NPt;
t = linspace(0,T,length(Po));
E = (1/dt)*(linspace(-pi,pi,length(Po)));
%%
Po = (1-cos(2*pi*t/T)).*Po;
Pe = fftshift(fft((Po/T)));
%%
figure(2);subplot(2,1,1);plot(t,real(Po));
title('Correlation Function ');xlabel('Time');
figure(2);subplot(2,1,2);plot(E,log(fliplr(abs(Pe))),'r');
title('Energy Spectrum');xlabel('Energy');ylabel('Power');
axis([-210 0 -17 5]);
pause(1);

%%-------------------------------------------------------------------------
%% Analytic method: For Even Solutions (Even Wave functions)
%%
z0 = b*sqrt(2*M*V0);
z = 0:0.01:20*pi;
y1 = tan(z);
y2 = sqrt((z0./z).^2 - 1);
figure(3);subplot(2,1,1);plot(z,y1,z,y2);
hold on;
plot(z,0*z,'r');
axis([0 45 0 35]);
title('tan(z)  =  [(z_0/z)^2 - 1]^{1/2}');
crss_n = [1.5 4.5 7.6 10.8 13.83 16.9 20.0 23.0 26.1 29.1 32.2 35.2 38.2 41.1]; 
%% ^-- get these values by looking at the graph (approx)
g =  inline('tan(z) - sqrt((z0/z).^2 - 1)','z','z0');
for nrn = 1:14
zn(nrn) = fzero(@(z) g(z,z0),crss_n(nrn));
end
figure(3);subplot(2,1,1);hold on;plot(zn,tan(zn),'rx');
q = zn/b;
Em = ((q.^2)/(2*M))-V0;
%%
for nrn = 1:length(Em),
    figure(3);subplot(2,1,2);hold on; 
    plot([Em(nrn),Em(nrn)],[-17,6]); 
end
%%
figure(3);subplot(2,1,2);
plot(E,log(fliplr(abs(Pe))),'r');hold on;
title('Energy Spectrum (Blue: Even solutions)');
xlabel('Energy');ylabel('Power');
axis([-210 0 -17 5]);
%%
%%-------------------------------------------------------------------------

