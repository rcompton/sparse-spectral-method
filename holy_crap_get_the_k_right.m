%
% Test out the different time propagators
%
clear all;close all;clc;
scp_so_finitewell_ORIGINAL;
Pe0 = Pe;
save scp_orig;


scp_so_finitewell;
Pe1 = Pe;
save scp_diffk;

%% red is original..
load('scp_orig','Pe0');
load('scp_diffk','Pe1');
figure()
plot(abs(Pe1));
hold on;
plot(abs(Pe0),'r--');

norm(Pe0 - Pe1)/norm(Pe0)