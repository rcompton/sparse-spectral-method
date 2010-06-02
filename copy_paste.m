%% Finite Square Well: Split-Operator Method
%% Ref: M. D. Feit, J. A. Fleck, Jr., and A. Steiger, 
%%      "Solution of the Schrodinger Equation by a Spectral Method",
%%      Journal of Computational Physics 47, 412-433 (1982).
%% Ref: D. J. Griffiths, 
%%      "Introduction to Quantum Mechanics",
%%      ISBN 0-13-124405-1
%% 
clear all; close all;clc;
a = 48;  %% Length
M = 1/2; %% Mass
N = 512;
x = linspace(-a/2,a/2,N); x = x'; %work with many rows and 1 column.
k = N*linspace(-1/2,1/2,N); k = k';
dt = 1e-3; %% Time step (dt < pi/(3*dVmax)...)
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
title('Potential');
plot(x,V,'r');hold on;plot(x,max(abs(real(V)))*abs(Phi0c));hold off; pause(1);
%%
%the method in it's pure form applys GK 2x
%but thats' the same and GK2 once in a loop...
GK = fftshift(exp(-(1i*dt/(4*M))*((2*pi/a)^2)*(k.^2))); %% dt/2 kinetic energy propagator
GK2 = fftshift(exp(-(1i*dt/(2*M))*((2*pi/a)^2)*(k.^2))); %% dt kinetic energy propagator
GV = exp(-1i*dt*V); %% Potential spatial interaction
% plot((-(dt/(4*M))*((2*pi/a)^2)*(k.^2)));
% plot(-dt*V);
%%
iPhi = fft(Phi0);
Phi = ifft(iPhi.*GK); %propagate dt/2 to get in line with future loop iters...
Phi = GV.*Phi;
NPt = 50000; %How did you pick? dEmin < pi/T is minimum seperation that can resolve.
Pt = zeros(1,NPt);
En = -105.99;  %% Energy eigenvalue. what the hell? how???
T = dt*NPt;
t = linspace(0,T,NPt);
%wt = (1-cos(2*pi*t/length(t))); %I don't think I'll need the Hanning window.
%uns = 0; 

%
% Propagate all the timesteps
% (lots of computational work...)
%
for nrn = 1:NPt
    %note, Phi was already propagated dt/2 before the loop.
    iPhi = fft(Phi);
    Pt(nrn) = trapz(x,Phi0c.*ifft(iPhi.*GK)); %propagate by dt/2 to line up with paper.
    Phi = ifft(iPhi.*GK2); %propagate by dt
%%
%     %%
%     %
%     % Compute the eigenfunctions
%     % (which I don't care about)
%     %
%     unl = Phi*wt(nrn)*exp(1i*En*nrn*dt);
%     if nrn > 1
%         una = (unp + unl)*dt/2; %% Trapezoidal area
%         uns = uns + una;  %% Explicit trapezoidal integration
%         if mod(nrn,1000) == 0
%             figure(1);
%             subplot(4,1,3);plot(x,abs(una)); 
%             title('\int_t^{t+dt} \Phi_x(t) w(t) exp(i*E_n*t) dt');
%             xlabel('x');axis tight;
%             subplot(4,1,4);plot(x,real(uns)); 
%             title('Eigen function @ E = -105.99');
%             xlabel('x');ylabel('Amp');axis tight;
%         end
%     end
%     unp = unl;
%     %%
%     if mod(nrn,1000) == 0
%         figure(1);
%         subplot(4,1,1);plot(x, real(Phi),'r'); 
%         title(['\Phi_x  t=',num2str(t(nrn))]); 
%         xlabel('x');ylabel('Amp');axis tight;
%         subplot(4,1,2);plot(k,fftshift(real(iPhi)),'r-');
%         title('\Phi_k'); 
%         xlabel('k');ylabel('k-space Amp');axis tight;
%         pause(0.2);
%     end
%%
%%
    Phi = GV.*Phi;    
end
iPhi = fft(Phi);
Phi = ifft(iPhi.*GK);%line up with paper's timestepping by dt/2 propagation
%%
estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
T = dt*NPt;
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

