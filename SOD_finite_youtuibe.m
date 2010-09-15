%SQUAREWELL.M

clear all; close all; clc;

%L = 2*pi; % Length of interval
n = 2*64; % number of divisions of x
dx = .32/3;
L = n*dx;
%x = 0:dx:n;
x = linspace(-L, L, n); 


u = zeros(1, n);
ut = zeros(1, n);
fps = 30; %frames per second
tmax = 15*fps;
dx = x(21) - x(20); 
steps = 1000; % steps calculated per frame dt < .25*dx^2
dt = 1 / steps; %time step

%Parameters and initial consditions
%u0, V, a
t = 0;
u = exp(-2*(x-L/3).^2 + 0*1i*x);

%V = zeros(1,n);
V = .2*(x>2*L/3) + .1*(x<0);

M = 50;

%normalization
pdf = u.*conj(u);
I = sum(pdf) * dx;
u = u / sqrt(I);

%animation capture setup
fig1 = figure(1);
axis([0, L, 0, 1]);
aviobj = avifile('swe1d3.avi', 'FPS', fps);
f = 0;

%a = .01

% u_t = -iHu = -i((-1/2M)*lap + V)u
H = inline( '(-1/(2*M))*4*del2(Phi,dx) + V.*Phi','M', 'V', 'dx', 'Phi');
% 
% k = n*linspace(-1/2,1/2,n); k = k';
% k = -n/2:n/2-1;
% k = k';
% k = k*2*pi/abs(max(x)-min(x));
% k = fftshift(k);
% H = inline( '(-1/(2*M))*ifft(-(k.^2).*fft(Phi)) + V.*Phi','M','V','k','Phi');


%Euler to get u1
unm1 = u;
un = u - 1i*dt*H(M,V,dx,unm1);

while t <= tmax
    %ut = -1i*((-1/(2*M))*4*del2(u,dx) +V.*u); % i*a * u_xx + i*V*u
    %ut = -1i*H(M,V,dx,u);
    %u = u + ut * dt; %Euler's method        
    u = unm1 - 2i*dt*H(M,V,dx,un);
    un = u;
    unm1 = un;
    
    % plot u, V, pdf
    f = mod(f+1, steps);
    if f==1
        plot(x,V,'color','g','linewidth',2);
        hold on
        plot(x,real(u),'--','color','r')
        plot(x,imag(u),'--','color','b')
        pdf = u.*conj(u);
        plot(x,pdf,'-','linewidth',1,'color','k','linewidth',2);
        axis([min(x), max(x), -1, 4]);
        %axis off
        hold off
        pause(.001)
        F = getframe(fig1);
        aviobj = addframe(aviobj,F);
    end
    t = t + dt;
end

close(fig1);
aviobj = close(aviobj);
