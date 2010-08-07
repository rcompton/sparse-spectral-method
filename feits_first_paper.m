%% Finite Square Well: Split-Operator Method
%% Ref: M. D. Feit, J. A. Fleck, Jr., and A. Steiger,
%%      "Solution of the Schrodinger Equation by a Spectral Method",
%%      Journal of Computational Physics 47, 412-433 (1982).
%% Ref: D. J. Griffiths,
%%      "Introduction to Quantum Mechanics",
%%      ISBN 0-13-124405-1
%%
clear all; close all; clc;
% 
M = 1/2; % particle mass
% N = 512; % number of space points
% x = linspace(-a/2,a/2,N); x = x'; % work with columns

% Define spatial grid vairables
dx = 0.0825;
N = 512;
L0 = dx*N;
xmax = L0/2;
xmin = -L0/2;
x = linspace(xmin,xmax,N); x = x';
k = N*linspace(-1/2,1/2,N); k = k'; % space frequency variables

% weird timescale stuff
Nt = 50*16384;
dt = 5.73;
T = dt*Nt;
T = 9.388e02;
dt = T/Nt;

t = linspace(0,T,Nt); % high resolution timescale
Pt = zeros(1,Nt); % autocorrlelation function has lots of zeros

%wt = (1-cos(2*pi*t/length(t))); %window of the high res scale


%% Potential

% Define potential
% For some weird physics reason it's mostly zeros...
% k0 = -132.7074997;
% k2 = 7;
% k3 = 0.5;
% k4 = 1;
% x1 = 3.813;
% x2 = -4.112;
% V = (k0 - k2*x.^2 + k3*x.^3 + k4*x.^4).*(x<x1).*(x2<x);

% % Morse Potential
De = 0.01688;
beta = 1.47612;
re = 3.28892;
%x = linspace(0.1,50,N);
V = De*(1 - exp(-beta.*(x - re)) ).^2 - De;
V = V.*(V < .05);
%V = V';

fprintf('dt %f , pi/min(V) %f \n', dt, pi/max(abs(V)) );

%% Define initial wavefunction
a = 1.9; % shift
sigmah = 0.87; % spread
Phi0 = exp(-((x-a).^2)/(2*sigmah^2)) + exp(-((x+a).^2)/(2*sigmah^2));

%Phi0 = exp(-((x-3.0)/.25).^2);

Phi0 = Phi0/norm(Phi0);
Phi0c = conj(Phi0); % real(psi0)- i*imag(psi0);

%% Make the propagators
%Note to self, the sprintf trick works to get scalars in scope
%but what about a vector?
%also note, GK as a function makes everything take forever!

%GK = inline(sprintf('fftshift(exp(-(i*dt/(4*%f))*((2*pi/%f)^2)*(k.^2)))',M,a)','dt','k'); %% dt/2 kinetic energy propagator
%GV = inline(sprintf('exp(-1i*dt*V)'),'dt','V'); %% Potential spatial interaction

GK = fftshift(exp(-(i*dt/(4*M))*((2*pi/a)^2)*(k.^2))); %% dt/2 kinetic energy propagator
GV = exp(-i*dt*V); %% Potential spatial interaction

%% Propagate all the timesteps
% ***One timestep per iteration!***
Phi = Phi0;
iPhi = fft(Phi);
%fprintf('we are going to make %d timesteps now:\n', max(size(dts)));
tic
for nrn = 1:Nt
    % momentum space propagation
%    iPhi = iPhi.*GK(dts(nrn),k);
    iPhi = iPhi.*GK;
    % move into physical space and apply potential operator
    Phi = ifft(iPhi);
%    Phi = GV(dts(nrn),V).*Phi;
    Phi = GV.*Phi;
    % move into momentum space and propagate again
    iPhi = fft(Phi);
    %iPhi = iPhi.*GK(dts(nrn),k);
    iPhi = iPhi.*GK;
    % move back into physical space and record P(t) at the sample point
    Phi = ifft(iPhi);
%    Pt(t_idxs(nrn)) = trapz(x, Phi0c.*Phi);
    Pt(nrn) = trapz(x,Phi0c.*Phi);
   
    if mod(nrn,5000) == 0
        fprintf('percent done %f\n', 100*nrn/Nt);
    end
end
toc
%%
estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
%E = (1/dt_small)*(linspace(-pi,pi,length(Pt)));
E = (1/dt)*linspace(-pi,pi,length(Pt));
Po = (1-cos(2*pi*t/T)).*Po;
%% Two ways...
% This way works for sure...
Pe = fft(Po);
Pe = fftshift(Pe)/T;

%%Below is the FPC_AS way,
%%it also works well but I don't
%%know what to do about the conjugation.
% 
% practice_points = sort(randsample(1:NPt, NPt/5.5));
% Pder = zeros(1,NPt);
% Pder(practice_points) = Po(practice_points);
% 
% % FPC_AS A_operator class
% A = A_operator( @(z) pifft(z, find(Pder)), @(z) pfft(z, find(Pder), NPt) );
% mu = 1e-10;
% [Pe, ~] = FPC_AS(NPt, A, nonzeros(Pder), mu);
% % Hurr durr scale by T
% Pe = Pe*sqrt(NPt)/T;
% Pe = conj(Pe);%??????????SHIT SHIT SHIT SHIT?????????????!!!!!!!!!!!
% %This...
% Pe = fftshift(Pe);


%% The exact according to Chen

% omegahe = .00632;
% 
% n_Es = 30;
% nPhalf = 1:n_Es + .5*ones(1,n_Es);
% En = omegahe*nPhalf - (nPhalf.^2)*(omegahe^2)/(4*De) - De;
% deltaE = (E(10)-E(9))
% Eeexact = zeros(size(E));
% for i=1:(n_Es)
%     Eeexact((En(i)+deltaE >= E) & (E >= En(i)-deltaE) )= 1;
% end

%% Plot

figure();
plot(E, log(fliplr(abs(Pe))));

figure();
plot(E,abs(Pe));

