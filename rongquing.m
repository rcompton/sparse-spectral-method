%
% Party script
%
%

clear all;
clc;

Nt = 2^13;
Nr = 2^12;

%.1 to 50au...
r = 0.1:1/Nr:50;

%
%Does is integrate to 1?
%
trapz(r, abs(psi_not(r)).^2)
