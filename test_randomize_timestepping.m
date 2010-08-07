%
% Script to test a randomized timestepping scheme
%
% We solve y' = f(y) at random locations and observe
% the errors
%
clear all;close all; clc;

N = 1421;
T = 17;
dt = T/N;

f = inline('y+y^-2', 'y');
y0 = 1;

%
% Euler's method with uniform steps
%
y = zeros(1,N);
y(1) = y0;
for i=1:N-1
    y(i+1) = y(i) + f(y(i))* dt;
end

figure();
plot(y, 'bd');

%
% Random steps for Euler's method
%
num_steps = N/2;
step_idxs = sort(randsample(1:N, num_steps));
step_idxs(1) = 1;

gap_sizes = [step_idxs 1e6] - [1e6 step_idxs];
gap_sizes = gap_sizes(2:num_steps);
gap_sizes = gap_sizes*dt;


yR = zeros(1,num_steps);

yR(1) = y0;
for i=1:num_steps-1
    yR(i+1) = yR(i) + f(yR(i))*gap_sizes(i);
end

yPlot = zeros(1,N);
yPlot(step_idxs) = yR;

hold on;
plot(yPlot,'rx');

fprintf('random error: %f \n', norm(y(step_idxs)-yR));

