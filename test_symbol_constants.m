clear all;close all;clc;


%% From Trefelthen's Book
% set up parameters

N = 24163;

%side length (we want an xgrid of -L,L)
L = 2.4*pi;

%dx = 2*pi/N; %relation: pi/dx = N/2, the x-space grid is [0,2*pi]
%x = dx*(1:N); %[0,2*pi] in N points
%x = L*(x-pi)/pi; %[-L,L] in N points

x = linspace(-L,L,N);
dx = median(diff(x));

fprintf('dx two ways %f %f \n', dx, L*2/N);
fprintf('xmin and xmax %f %f \n', min(x) ,max(x));

%k = [0:N/2-1 0 -N/2+1:-1]; %???
k = (pi/L)*(-N/2:N/2-1); % pi/L is dk
k = fftshift(k);

%% So... if it's periodic then you know what to do.
%if it's not, then you can make it periodic by padding with zeros?
periodic_trash = inline('exp(sin(x))-exp(cos(x).^2)','x');
v = periodic_trash(x);

%physical space differentiation
vprime = diff(v)/dx;

%k-space multiply
v_hat = fft(v);
w_hat = 1i.*k.*v_hat;
w = real(ifft(w_hat));

close all;
plot(x,w,'rx');hold on;plot(x(2:length(x)),vprime);


%% again...

%This is periodic as it has compact support
%and is zero at endpoints of intervals
v = exp(-x.^2);

%physical space differentiation
vprime = diff(v)/dx;
vprimeprime = diff(vprime)/dx;

%k-space multiply
v_hat = fft(v);
w_hat = ((1i.*k).^2).*v_hat;
w = real(ifft(w_hat));

close all;
plot(x,w,'rx');hold on;plot(x, 4*del2(v,dx));

%% and for 2d...

N = 200;
[xx yy] = meshgrid( linspace(-L,L,N), linspace(-L,L,N) );
dx = median(xx);
dy = dx;

[kk ll] = meshgrid( -N/2:N/2-1, -N/2:N/2-1); % pi/L is dk
kk = (pi/L)*fftshift(kk);
ll = (pi/L)*fftshift(ll);

v2d = exp(-xx.^2 - yy.^2);

%k-space multiply
v2d_hat = fft2(v2d);
w2d_hat = ( -kk.^2 - ll.^2 ).*v2d_hat;
w2d = real(ifft2(w2d_hat));


figure();
plot3(xx,yy, 4*del2(v2d,dx,dy),'r--');
hold on;
plot3(xx,yy, w2d);

figure()
plot3(xx,yy, 4*del2(v2d,dx,dy) - w2d);
