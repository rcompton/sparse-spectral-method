%SQUAREWELL.M
%2D

clear all; close all; clc;

%L = 2*pi; % Length of interval
n = 128; % number of divisions of x
L = 32.5;

[xx yy] = meshgrid(linspace(-L/2, L/2, n), linspace(-L/2, L/2, n));

dx = median(diff(linspace(-L/2,L/2,n)));
dy = dx;

fps = 30; %frames per second
%tmax = 5*fps;
tmax = 314;
dt = 0.05 %time step, unstable at 1/200 with SOD!

%each frame is 1 unit of time.
%for a chebyshev video you'll want 1 (or less) steps per frame
steps = 1.0/dt; % steps calculated per frame dt < .25*dx^2 ?Depends on M

M = 15;

%2D well
%Uh, I'm kinda dumb and did this a hard way
depth = 1.8;
width = L/4.5; %density of states is inverse to this...
V = depth*(xx<-width) + depth*(yy<-width) - depth*(yy<-width).*(xx<-width);
V = V + rot90(V,2) -depth*(xx>width).*(yy<-width) -depth*(xx<-width).*(yy>width);

%Henion Hines
% lambdah = .118034;
% V = (M/22)*(xx.^2 + yy.^2) + lambdah*xx.*(yy.^2 - 0.33*xx.^2);
% V = V.*(V < 2.5) + 2.5*(V>2.5);
% 
% %2D double well
% k0 = -132.7074997;
% k2 = 7;
% k3 = 0.5;
% k4 = 1;
% x1 = 3.813;
% x2 = -4.112;
% V = k0 - k2*(xx.^2 + yy.^2) + k3*(xx.^3 + yy.^3) + k4*(xx.^4+yy.^4);
% V = V.*(xx.^2 + yy.^2 < 4.1);
% 
% V = V'/10 + min(V(:));


%plot3(xx,yy,V,'r.')
%surf(xx,yy,V)



%initial psi
%psi0 = exp(-((xx+0.85).^2 + (yy-1.0).^2)) + exp(-(((xx-0.85)/0.5).^2 +((yy+1.0)/0.5).^2));
%psi0 = exp(-((xx/0.5+1.5).^2 + (yy/0.5 -1.5).^2)) + exp(-((xx/0.3-2).^2 + (yy/0.3 +2.8).^2));
a = L/16; % shift
sigmah = .7; % spread
psi0 = exp(-((xx+a).^2 + (yy-a).^2)/(2*sigmah^2)) + exp(-((xx-a).^2 + (yy+a).^2 )/(2*sigmah^2));

%normalize
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
aviobj.Quality = 30;

%H, dx=dy
H = inline( '(-1/(2*M))*4*del2(Phi,dx,dx) + V.*Phi','M', 'V', 'dx', 'Phi');

%spectral style differencing.
% k = -n/2:(n/2-1);
% k = k.*pi./L;
% k = fftshift(k);
[kk ll] = meshgrid( -n/2:(n/2-1), -n/2:(n/2-1)); % pi/L is dk
kk = (2*pi/L)*fftshift(kk);
ll = (2*pi/L)*fftshift(ll);

% Differencing functions
% for normalized H you need to estimate minE etc.
%Hspec = inline( '(-1/(2*M))*ifft2(-(kk.^2 +ll.^2).*fft2(Phi)) + V.*Phi','M','V','kk','ll','Phi');
Hspec = inline( '(-1/(2*M))*specdiff2d(Phi,xx,yy) + V.*Phi','M','V','xx','yy','Phi');
dE = (pi^2)/(2*M*dx^2) + max(V) - min(V)
maxE = (pi^2)/(2*M*dx^2) + max(V)
minE = min(V)
Hnormspec = inline( '(2/dE)*((-1/(2*M))* ifft2(-(kk.^2 +ll.^2).*fft2(Phi)) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','kk','ll','Phi');
%Hnormspec = inline( '(2/dE)*((-1/(2*M))* specdiff2d(Phi,xx,yy) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','xx','yy','Phi');


%check it, the kk and ll
%close all;
%surf(xx,yy,V,'EdgeColor','none');
%plot3(xx,yy,abs(Hspec(M,V,xx,yy,psi0))-abs(H(M,V,dx,psi0)),'r-');
%%


dE = (pi^2)/(2*M*dx^2) + max(V(:)) - min(V(:))
maxE = (pi^2)/(2*M*dx^2) + max(V(:))
minE = min(V(:))


%
NPt = ceil(tmax/dt);

%set up a randomized time grid
num_samples = ceil(NPt/1);
stream = RandStream('mrg32k3a');
sample_points = sort(unique([1 randsample(stream,1:NPt, num_samples)]));
num_samples = length(sample_points);

ts = sample_points*dt;
dts = diff(ts);

Pt = zeros(1,ceil(tmax/dt));

%Euler to get u1
%unm1 = psi0;
%un = psi0 - 1i*dt*H(M,V,dx,unm1);
%un = u - 1i*dt*Hnormspec(M,V,dE,minE,kk,ll,unm1);


% Store them all.
%psis = zeros(n,n,NPt);
%psis(:,:,1) = psi0;
psisnrnm1 = psi0;

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
    Tkpsis(:,:,1) = psisnrnm1;
    Tkpsis(:,:,2) = Hnormspec(M,V,dE,minE,kk,ll,psisnrnm1);
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
    
    psisnrnm1 = chebsum;
    
    % Plotting
    if mod(nrn,22)==0
        plot3(xx,yy,V,'r:');
        hold on;
        surf(xx,yy,abs(chebsum),'EdgeColor','none');        
        axis([min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:)) -0.05 0.45]);

        %pause(.01);
        hold off;
        F = getframe(figgn);
        aviobj = addframe(aviobj,F);
        tmax
        t
        
        %plot to file
        figure(2); set(2,'visible', 'off');
        plot3(xx,yy,V,'r:');
        hold on;
        surf(xx,yy,abs(chebsum),'EdgeColor','none');        
        axis([min(xx(:)) max(xx(:)) min(yy(:)) max(yy(:)) -0.05 0.45]);
        hgexport(2,strcat('graph',num2str(nrn)) );close(2);
    end
    
    t = t+dt;
    norm(chebsum);
    derp = chebsum.*conj(psi0);
    Pt(nrn) = trapz(trapz(derp))*dx*dy;
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

