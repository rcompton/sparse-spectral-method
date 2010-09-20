%SQUAREWELL.M

clear all; close all; clc;

%L = 2*pi; % Length of interval
n = 512; % number of divisions of x
L = 50;

%x = 0:dx:n;
x = linspace(0, L, n); 
dx = x(21) - x(20)


u = zeros(1, n);
ut = zeros(1, n);
fps = 30; %frames per second
tmax = 5*fps;
steps = 200; % steps calculated per frame dt < .25*dx^2 ?Depends on M
dt = 100 / steps; %time step



%Parameters and initial consditions
%u0, V, a
t = 0;
%u = exp(-2*(x-L/3).^2 + 0*1i*x);

%V = zeros(1,n);
%V = .02*(x+L/6).*(x>-L/6) + .1*(x<-L/6);
%V = .3*(x<-L/3) + .5*(x>2*L/3);
%V(1:10) = 10;
%V(end-10:end) = 10;

re = 3.2889;
De = 0.01688;
betah = 1.47612;
u = exp(-((x-3)/.25).^2);
V = De*(1-exp(-betah*(x-re))).^2 - De;
V = V.*(V < 0.5) + (V>0.5);



M = 50;
fprintf('stability condition %f \n', (1/(2*M))*dt/dx^2);
dt
.25*dx^2

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

%spectral style differencing.
%k = n*linspace(-1/2,1/2,n); k = k';

k = -n/2:n/2-1;
%k = k';
%k = k*1*pi/abs(max(x)-min(x));
k = k*1*pi/L;
k = fftshift(k);
Hspec = inline( '(-1/(2*M))*ifft(-(k.^2).*fft(Phi)) + V.*Phi','M','V','k','Phi');


u0 = u;
Pt = zeros(1,ceil(tmax/dt));

%Euler to get u1
unm1 = u;
un = u - 1i*dt*H(M,V,dx,unm1);

counter = 1;
while t <= tmax
    %ut = -1i*((-1/(2*M))*4*del2(u,dx) +V.*u); % i*a * u_xx + i*V*u
    %ut = -1i*H(M,V,dx,u);
    %u = u + ut * dt; %Euler's method        
    
    %SOD without inline function.
    %Hu = (-1/(2*M))*([0 0 u] - 2*[0 u 0] + [u 0 0])/(dx^2);
    %Hu = Hu(2:end-1) + V.*u;
    %u = unm1 - 2i*dt*Hu;
    
    
    %u = unm1 - 2i*dt*H(M,V,dx,un); %SOD
    %u = unm1 - 2i*dt*Hspec(M,V,dx,un); %SOD Fourier

    %un = u;
    %unm1 = un;
    
    %cross fingers
    if mod(steps,700)==0
        u = chebystep(M,dx,V,dt,u,true);
    else
        u = chebystep(M,dx,V,dt,u,false);
    end
    %u = chebyshev_apply(dx, M, V, k, dt, u);
    
    
    % plot u, V, pdf
    f = mod(f+1, steps/5);
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
        pause(.0001)
        F = getframe(fig1);
        aviobj = addframe(aviobj,F);
    end
    t = t + dt;

    
    Pt(counter) = trapz(x,u.*conj(u0));
    counter = counter + 1;
end

close(fig1);
aviobj = close(aviobj);

% back to old shit
ts = linspace(0,tmax,length(Pt));
wt = (1-cos(2*pi*ts/length(ts)));

estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
E = (1/dt)*(linspace(-pi,pi,length(Pt)));
Po = (1-cos(2*pi*ts/tmax)).*Po;

Pe = fft(Po);
Pe = fftshift(Pe)/tmax;

figure();
plot(Pe);

