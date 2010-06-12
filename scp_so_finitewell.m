%% Finite Square Well: Split-Operator Method
%% Ref: M. D. Feit, J. A. Fleck, Jr., and A. Steiger, 
%%      "Solution of the Schrodinger Equation by a Spectral Method",
%%      Journal of Computational Physics 47, 412-433 (1982).
%% Ref: D. J. Griffiths, 
%%      "Introduction to Quantum Mechanics",
%%      ISBN 0-13-124405-1
%% 
clear all; close all; clc;

a = 48;  % space size
M = 1/2; % particle mass
N = 512; % number of space points
x = linspace(-a/2,a/2,N); x = x'; % work with columns
k = N*linspace(-1/2,1/2,N); k = k'; % frequency vairables
NPt = 50000 % number of time steps
dt = 1e-3 % smallest possible time step
Pt = zeros(1,NPt); % autocorrlelation function
T = dt*NPt % max time
t = linspace(0,T,NPt); %high resolution timescale
wt = (1-cos(2*pi*t/length(t))); %window of the high res scale

%% randomized timestepping
p = 1
num_time_samples = NPt*p;
t_samp = sort(randsample(t, num_time_samples));
dts = [t_samp 0] - [0 t_samp];
dts = dts(1:num_time_samples-1); %careful here...

%% Potential
V0 = 200;
V = zeros(length(x),1) - V0; % 1*((2*x).^2 - (0.6*a)^2);% iPhi = fft(Phi0);
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
%figure(1);set(gcf,'position',[37 208 538 732]);
%plot(x,V,'r');hold on;plot(x,max(abs(real(V)))*abs(Phi0c));hold off; pause(1);

%% Make the propagators
%Note to self, the sprintf trick works to get scalars in scope
%but what about a vector?
%also note, GK as a function makes everything take forever!
GK = inline(sprintf('fftshift(exp(-(i*dt/(4*%f))*((2*pi/%f)^2)*(k.^2)))',M,a)','dt','k'); %% dt/2 kinetic energy propagator
GV = exp(-1i*dt*V); %% Potential spatial interaction

% plot((-(dt/(4*M))*((2*pi/a)^2)*(k.^2)));
% plot(-dt*V);
%%

% Propagate all the timesteps
% ***One timestep per iteration!***
Phi = Phi0;
iPhi = fft(Phi);
for nrn = 1:NPt
    % momentum space propagation
    iPhi = iPhi.*GK(dts(nrn),k);
    % move into physical space and apply potential operator
    Phi = ifft(iPhi);
    Phi = GV.*Phi;
    % move into momentum space and propagate again
    iPhi = fft(Phi);
    iPhi = iPhi.*GK(dts(nrn),k);
    % move back into physical space and record P(t)
    Phi = ifft(iPhi);
    Pt(nrn) = trapz(x, Phi0c.*Phi);
end


%%
estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
T = dt*NPt;
t = linspace(0,T,length(Po));
E = (1/dt)*(linspace(-pi,pi,length(Po)));
Po = (1-cos(2*pi*t/T)).*Po;
%% Two ways...
%Pe = fft(Po);
Pe = fft_via_fpca(Po, randsample(1:NPt, NPt/3));

Pe = fftshift(Pe)/T;
%%

% figure(2);subplot(2,1,1);plot(t,real(Po));
% title('Correlation Function ');xlabel('Time');
% figure(2);subplot(2,1,2);plot(E,log(fliplr(abs(Pe))),'r');
% title('Energy Spectrum');xlabel('Energy');ylabel('Power');
% axis([-210 0 -17 5]);
% pause(1);

%%-------------------------------------------------------------------------
%% Analytic method: For Even Solutions (Even Wave functions)
%%
z0 = b*sqrt(2*M*V0);
z = 0:0.01:20*pi;
y1 = tan(z);
y2 = sqrt((z0./z).^2 - 1);
%figure(3);subplot(2,1,1);plot(z,y1,z,y2);
%hold on;
%plot(z,0*z,'r');
%axis([0 45 0 35]);
%title('tan(z)  =  [(z_0/z)^2 - 1]^{1/2}');
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

