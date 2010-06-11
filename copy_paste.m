%% Split-Operator Method
%% Ref: M. D. Feit, J. A. Fleck, Jr., and A. Steiger, 
%%      "Solution of the Schrodinger Equation by a Spectral Method",
%%      Journal of Computational Physics 47, 412-433 (1982).



%% 
clear all;close all;clc;

% Define spatial grid vairables
dx = 0.825
Nx = 512
L0 = dx*Nx
xmax = L0/2
xmin = -L0/2
x = linspace(xmin,xmax,Nx);

% Define potential
% For some weird physics reason it's mostly zeros...
k0 = -132.7074997;
k2 = 7;
k3 = 0.5;
k4 = 1;
x1 = 3.813;
x2 = -4.112;
V = (k0 - k2*x.^2 + k3*x.^3 + k4*x.^4).*(x<x1).*(x2<x);

% Define uniform time domain
% WTF?!
dt = .00573 %% Time step (dt < pi/(3*dVmax)...)
Nt = 50384

% Define initial wavefunction
a = 1.9; % shift
sigmah = 0.87; % spread
psi0 = exp(-((x-a).^2)/(2*sigmah^2)) + exp(-((x+a).^2)/(2*sigmah^2));
psi0c = conj(psi0); % real(psi0)- i*imag(psi0);

% Mass
M = 1/2

% Define spectral vairables
k = Nx*linspace(-1/2,1/2,Nx);

%% Define the spectral operators.
% the method in it's pure form applys GK twice
% but thats' the same as an offset and GK2 once in a loop...
% why do we shift?
% I don't know.
GK = fftshift(exp(-(1i*dt/(4*M))*((2*pi/L0)^2)*(k.^2))); %% dt/2 kinetic energy propagator
GV = exp(-1i*dt*V); %% Potential spatial interaction


%% The major computational loop
ipsi = fft(psi0);
psi = ifft(ipsi.*GK); %propagate dt/2 to get in line with future loop iters...
psi = GV.*psi;

Pt = zeros(1,Nt);
T = dt*Nt;
t = linspace(0,T,Nt);

%
% Propagate all the timesteps
% ***One timestep per iteration!***
%
psi = psi0;
ipsi = fft(psi);
for nrn = 1:Nt
    % momentum space propagation
    ipsi = ipsi.*GK;
    % move into physical space and apply potential operator
    psi = ifft(ipsi);
    psi = GV.*psi;
    % move into momentum space and propagate again
    ipsi = fft(psi);
    ipsi = ipsi.*GK;
    
    % move back into physical space and record P(t)
    psi = ifft(ipsi);
    Pt(nrn) = trapz(x, psi0c.*psi);
    
    
end



%%
estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
T = dt*Nt;
t = linspace(0,T,length(Po));
E = (1/dt)*(linspace(-pi,pi,length(Po)));
%% Draw the energy spectrum with classical FFT methods
Pe = fftshift(fft(((1-cos(2*pi*t/T)).*Po/T)));
Peww = fftshift( fft(Po/T) ); %no windowing in compressed sensing?

% figure(4);
% title('P(E) with windows, blue. Without windows, red');
% plot(E, Pe, 'b');
% hold on
% plot(E, Peww, 'r');

%% Compressed sensing reconstruction of energy spectrum
% we are given the whole Pt.
% Here we downsample it and then L1 reconstruct Pe
%
%
% 
% stream = RandStream('mrg32k3a');
% 
% %downsample the complete dataset
% num_samples = round(NPt/3.1234);
% sample_points = randsample(stream, 1:NPt, num_samples);
% 
% %Poder is the downsampled Pt
% %sparse freqs
% Poder = sqrt(NPt)*Pt(sample_points);
% 
% 
% %
% % We need to take an fft to get Pe
% %
% %A = A_operator( @(z) (sqrt(NPt))*R*ifft(z), @(z) (1/sqrt(NPt))*fft(R'*z));
% A = A_operator( @(z) pifft(z,sample_points), @(z) pfft(z,sample_points,NPt) );
% 
% % tiny mu corresponds to heavy weight on the fidelity term
% mu = 1e-10;
% 
% 
% % Call Wotao's code.
% [Peder, Out] = FPC_AS(NPt, A, Poder, mu);
% 
% %Peder = l1eq_pd(x0, A, At, Poder, 1e-3);
% 
% 
% hold on
% plot(E, Peder, 'g');

%%
% figure(2);subplot(2,1,1);plot(t,real(Po));
% title('Correlation Function ');xlabel('Time');
% figure(2);subplot(2,1,2);plot(E,log(fliplr(abs(Pe))),'r');
% title('Energy Spectrum');xlabel('Energy');ylabel('Power');
% axis([-210 0 -17 5]);
% pause(1);
% 
% %%-------------------------------------------------------------------------
% %% Analytic method: For Even Solutions (Even Wave functions)
% %%
% z0 = b*sqrt(2*M*V0);
% z = 0:0.01:20*pi;
% y1 = tan(z);
% y2 = sqrt((z0./z).^2 - 1);
% figure(3);subplot(2,1,1);plot(z,y1,z,y2);
% hold on;
% plot(z,0*z,'r');
% axis([0 45 0 35]);
% title('tan(z)  =  [(z_0/z)^2 - 1]^{1/2}');
% crss_n = [1.5 4.5 7.6 10.8 13.83 16.9 20.0 23.0 26.1 29.1 32.2 35.2 38.2 41.1]; 
% %% ^-- get these values by looking at the graph (approx)
% g =  inline('tan(z) - sqrt((z0/z).^2 - 1)','z','z0');
% for nrn = 1:14
% zn(nrn) = fzero(@(z) g(z,z0),crss_n(nrn));
% end
% figure(3);subplot(2,1,1);hold on;plot(zn,tan(zn),'rx');
% q = zn/b;
% Em = ((q.^2)/(2*M))-V0;
% %%
% for nrn = 1:length(Em),
%     figure(3);subplot(2,1,2);hold on; 
%     plot([Em(nrn),Em(nrn)],[-17,6]); 
% end
% %%
% figure(3);subplot(2,1,2);
% plot(E,log(fliplr(abs(Pe))),'r');hold on;
% title('Energy Spectrum (Blue: Even solutions)');
% xlabel('Energy');ylabel('Power');
% axis([-210 0 -17 5]);
% %%
% %%-------------------------------------------------------------------------

