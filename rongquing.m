%
% Party script
%
%

clear all;
clc;

Nt = 2^13;
dt = 4;

Nr = 2^12;

%.1 to 50au...
r = 0.1:1/Nr:50;

%
%Does is integrate to 1?
%
trapz(r, abs(psi_not(r)).^2)

%
% Advance some timesteps
%
clear all; close all;
Lnot = 48;  %% Length
M = 1/2; %% Mass
N = 512;
r = linspace(-Lnot/2,Lnot/2,N); r = r';
k = N*linspace(-1/2,1/2,N); k = k';
dt = 1e-3; %% Time step

V0 = 200;
V = zeros(length(r),1) - V0; % 1*((2*x).^2 - (0.6*a)^2);  %  
b = Lnot/16;
V(r<-b) = 0;
V(r>+b) = 0;

%figure()
%plot(r,V);


Psi0 = exp(-(5*(r-0*Lnot/128)).^2);

figure()
plot(r,Psi0)

psit = adv_one_step(r, Psi0, V, dt, M, Lnot);

for nt=1:3000
    psit = adv_one_step(r, psit, V, dt, M, Lnot);
end


hold on
plot(r,psit,'r')

