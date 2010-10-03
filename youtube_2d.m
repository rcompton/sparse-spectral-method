%SQUAREWELL.M
%2D

clear all; close all; clc;

%L = 2*pi; % Length of interval
n = ceil(sqrt(2102)); % number of divisions of x
L = 10;

[xx yy] = meshgrid(linspace(-L/2, L/2, n), linspace(-L/2, L/2, n));

dx = median(diff(linspace(-L/2,L/2,n)));
dy = dx;

fps = 30; %frames per second
%tmax = 5*fps;
tmax = 35;
dt = 1e-3 %time step, unstable at 1/200 with SOD!

%each frame is 1 unit of time.
%for a chebyshev video you'll want 1 (or less) steps per frame
steps = 1.0/dt; % steps calculated per frame dt < .25*dx^2 ?Depends on M

%2D well
%Uh, I'm kinda dumb and did this a hard way
depth = 1.2;
width = L/4;
V = depth*(xx<-width) + depth*(yy<-width) - depth*(yy<-width).*(xx<-width);
V = V + flipud(fliplr(V)) -depth*(xx>width).*(yy<-width) -depth*(xx<-width).*(yy>width);

%plot3(xx,yy,V,'r.')
%surf(xx,yy,V)

%normalization
psi0 = exp(-((xx+0.85).^2 + (yy-1.0).^2)) + exp(-(((xx-0.85)/0.5).^2 +((yy+1.0)/0.5).^2));
pdf = psi0.*conj(psi0);
I = sum(pdf(:)) * dx*dy;
psi0 = psi0 / sqrt(I);

%hold on;
%plot3(xx,yy,psi0,'b.')        set(1,'Visible','off')
%surf(xx,yy,psi0)


%animation capture setup
figgn = figure('Visible','off');
%set(h1,'Visible','off');
%axis([0, L, 0, 1]);
aviobj = avifile('swe1d3.avi', 'FPS', fps);
%f = 0;

%a = .01

M = 15;
%H, dx=dy
H = inline( '(-1/(2*M))*4*del2(Phi,dx,dx) + V.*Phi','M', 'V', 'dx', 'Phi');

%spectral style differencing.
k = -n/2:(n/2-1);
k = k.*pi./L;
k = fftshift(k);

[kk ll] = meshgrid( -n/2:(n/2-1), -n/2:(n/2-1)); % pi/L is dk
kk = (2*pi/L)*fftshift(kk);
ll = (2*pi/L)*fftshift(ll);

%
Hspec = inline( '(-1/(2*M))*ifft2(-(kk.^2 +ll.^2).*fft2(Phi)) + V.*Phi','M','V','kk','ll','Phi');

%check it, the kk and ll
%close all;
%plot3(xx,yy,abs(Hspec(M,V,kk,ll,psi0))-abs(H(M,V,dx,psi0)),'r-');



dE = (pi^2)/(2*M*dx^2) + max(V(:)) - min(V(:))
maxE = (pi^2)/(2*M*dx^2) + max(V(:))
minE = min(V(:))


Pt = zeros(1,ceil(tmax/dt));

%Euler to get u1
unm1 = psi0;
un = psi0 - 1i*dt*H(M,V,dx,unm1);
%un = u - 1i*dt*Hnormspec(M,V,dE,minE,k,dx,unm1);


counter = 1;
t=0;
while t <= tmax
    %ut = -1i*((-1/(2*M))*4*del2(u,dx) +V.*u); % i*a * u_xx + i*V*u
    %ut = -1i*H(M,V,dx,u);
    %u = u + ut * dt; %Euler's method
    
    %u = unm1 - 2i*dt*H(M,V,dx,un); %SOD
    %u = unm1 - 2i*dt*Hspec(M,V,dx,un); %SOD Fourier
    %
    %SOD without inline function.
    %Hu = (-1/(2*M))*([0 0 u] - 2*[0 u 0] + [u 0 0])/(dx^2);
    %Hu = Hu(2:end-1) + V.*u;
    
    %Fourier method (better)
    %Hu = (-1/(2*M))*ifft(-(k.^2).*fft(u));
    %Hu = Hu + V.*u;
    
    %Hu = H(M,V,dx,un);
    Hu = Hspec(M,V,kk,ll,un);
    
    u = unm1 - 2i*dt*Hu;
    un = u;
    unm1 = un;
    
    if mod(counter,33)==0
        plot3(xx,yy,V,'r.');
        hold on;
        surf(xx,yy,abs(u));
        %pause(.01);
        hold off;
        F = getframe(figgn);
        aviobj = addframe(aviobj,F);
    end
    
    %cross fingers
    %u = chebystep(M,dx,L,V,dt,u,true);
    %      if mod(steps,700)==0
    %          u = chebystep(M,dx,V,dt,u,true);
    %      else
    %          u = chebystep(M,dx,V,dt,u,false);
    %      end
    %
    
    % plot u, V, pdf
    %     f = mod(f+1, steps/5);
    %     if f==1
    %         plot(x,V,'color','g','linewidth',2);
    %         hold on
    %         plot(x,real(u),'--','color','r')
    %         plot(x,imag(u),'--','color','b')
    %         pdf = u.*conj(u);
    %         plot(x,pdf,'-','linewidth',1,'color','k','linewidth',2);
    %         axis([min(x), max(x), -1, 4]);
    %         %axis off
    %         hold off
    %         pause(.0001)
    %         F = getframe(fig1);
    %         aviobj = addframe(aviobj,F);
    %     end
    
    t = t + dt
    norm(u);
    derp = u.*conj(psi0);
    Pt(counter) = sum(derp(:))*dx*dy;
    counter = counter + 1;
end

close(figgn);
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
plot(abs(Pe));

