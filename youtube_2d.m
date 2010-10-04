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

% Differencing functions
% for normalized H you need to estimate minE etc.
Hspec = inline( '(-1/(2*M))*ifft2(-(kk.^2 +ll.^2).*fft2(Phi)) + V.*Phi','M','V','kk','ll','Phi');
dE = (pi^2)/(2*M*dx^2) + max(V) - min(V)
maxE = (pi^2)/(2*M*dx^2) + max(V)
minE = min(V)
Hnormspec = inline( '(2/dE)*((-1/(2*M))* ifft2(-(kk.^2 +ll.^2).*fft2(Phi)) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','kk','ll','Phi');


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
%un = u - 1i*dt*Hnormspec(M,V,dE,minE,kk,ll,unm1);


% Store them all.
NPt = ceil(tmax/dt);
psis = zeros(n,n,NPt);
psis(:,:,1) = psi0;

t=0;
for nrn=2:NPt
    dts(nrn) = dt;
    % SOD Fourier method.
%     u = unm1 - 2i*dt*Hspec(M,V,kk,ll,un);
%     un = u;
%     unm1 = un;
    
    % Do a chebystep
    %make maxk Tkpsis
    Tkpsis = zeros(n,n,100);
    Tkpsis(:,:,1) = psis(:,:,nrn-1);
    Tkpsis(:,:,2) = Hnormspec(M,V,dE,minE,kk,ll,psis(:,:,nrn-1));
    kunt = 2; nexttenjays = 1:.3:1.5;
    while max( abs(besselj(kunt,dE*dts(nrn-1)*nexttenjays)) ) > 1e-6
        kunt = kunt+1;
        %Tkpsis(:,k) = 2*Hnorm(M,V,dE,minE,dx,Tkpsis(:,k-1)) - Tkpsis(:,k-2);
        Tkpsis(:,:,kunt) = 2*Hnormspec(M,V,dE,minE,kk,ll,Tkpsis(:,:,kunt-1)) - Tkpsis(:,:,kunt-2);
    end
    
    %make the sum
    chebsum = besselj(0,dts(nrn-1)*dE)*Tkpsis(:,:,1);
    chebsum = chebsum + -2i*besselj(1,dts(nrn-1)*dE)*Tkpsis(:,:,2);
    for chebkunt=2:kunt-1
        chebsum = chebsum + 2*((-1i)^chebkunt)*besselj(chebkunt,dts(nrn-1)*dE)*Tkpsis(:,:,chebkunt+1);
    end
    
    psis(:,:,nrn) = chebsum;
    
    % Plotting
    if mod(nrn,33)==0
        plot3(xx,yy,V,'r.');
        hold on;
        surf(xx,yy,abs(chebsum));
        %pause(.01);
        hold off;
        F = getframe(figgn);
        aviobj = addframe(aviobj,F);
        tmax
        t
    end
    
    t = t+dt;
    %norm(chebsum);
    derp = chebsum.*conj(psi0);
    Pt(nrn) = sum(derp(:))*dx*dy;
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

