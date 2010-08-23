%
% Test out the different time propagators
%
clear all;close all;clc;
scp_so_finitewell_ORIGINAL;
Pe0 = Pe;
Pt0 = Pt;
save scp_orig;

clear all;
scp_so_finitewell;
Pe1 = Pe;
Pt1 = Pt;
save scp_diffk;

%% red is original..
load('scp_orig','Pe0','Pt0');
load('scp_diffk','Pe1','Pt1');
figure()
plot(abs(Pe1));
hold on;
plot(abs(Pe0),'r--');

norm(Pe0 - Pe1)/norm(Pe0)

% Note, the original code records the 1st Pt
% at time dt. The modified code records at time 0
% The original code gets to final time T+dt
% and records that as well. This is easily fixable
% and is not the problem...
%%
figure()
plot(real(Pt1));
hold on;
plot(real(Pt0),'r--');