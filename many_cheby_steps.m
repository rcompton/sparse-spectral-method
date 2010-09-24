%
% Propagate using Chebyshev.
%

close all;clear all;clc;
% Apply e^(lap + V) over a timestep dt

n = 512; % number of divisions of x (should be even!!)
L = 48;

x = linspace(-L/2, L/2, n);
dx = median(diff(x))


dt = .5;
tmax = 1700;
Nt = tmax/dt;

% re = 3.2889;
% De = 0.01688;
% betah = 1.47612;

%u = exp(-((x-3)/.25).^2);

%Morse
% V = De*(1-exp(-betah*(x-re))).^2 - De;
% V = V.*(V < 0.5) + (V>0.5);
% V = V';


%V = zeros(size(x));
%V = .7*(x<-L/3) + .8*(x>L/3);
%V = (2*x).^2;
V = 3.8*(x< -L/16) + 3.8*(x > L/16);
V = V';


% Wtf?
M = 50;

%animation capture setup
fig1 = figure(1);
fps = 30;
axis([0, L, 0, 1]);
aviobj = avifile('swe1d3.avi', 'FPS', fps);

dE = (pi^2)/(2*M*dx^2) + max(V) - min(V)
maxE = (pi^2)/(2*M*dx^2) + max(V)
minE = min(V)

eLL = max(x) - min(x);
k = -n/2:(n/2-1);
k = k*2*pi/(eLL);
k = fftshift(k)';


%psi0 = u;
%psi0 = exp(-((x-3)/.25).^2);
psi0 = exp(-(5*x).^2);
psi0 = psi0';
I = trapz(x,psi0.*conj(psi0));
psi0 = psi0 / sqrt(I);


Hnormspec = inline( '(2/dE)*((-1/(2*M))* specdiff(Phi,x) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','k','dx','x','Phi');

%figure maxk for a given timestep
maxk = 2;
nexttenjays = 1:.3:1.5;
while max( abs(besselj(maxk,dE*dt*nexttenjays)) ) > 1e-6 %&& abs(besselj(maxk,dt*(jdt+pi))) >1e-6
    maxk = maxk+1;
end

psis = zeros(n,Nt);
psis(:,1) = psi0;

for nrn = 2:Nt
      
    %% Do a chebystep
    %make maxk Tkpsis
    Tkpsis = zeros(n,100);
    Tkpsis(:,1) = psis(:,nrn-1);
    %Tkpsis(:,2) = Hnorm(M,V,dE,minE,dx,psi0);
    Tkpsis(:,2) = Hnormspec(M,V,dE,minE,k,dx,x,psis(:,nrn-1));
    kunt = 2;
    while max( abs(besselj(kunt,dE*dt*nexttenjays)) ) > 1e-6
        kunt = kunt+1;
        %Tkpsis(:,k) = 2*Hnorm(M,V,dE,minE,dx,Tkpsis(:,k-1)) - Tkpsis(:,k-2);
        Tkpsis(:,kunt) = 2*Hnormspec(M,V,dE,minE,k,dx,x,Tkpsis(:,kunt-1)) - Tkpsis(:,kunt-2);
    end
    kunt
    
    %make the sum
    chebsum = besselj(0,dt*dE)*Tkpsis(:,1);
    chebsum = chebsum + -2i*besselj(1,dt*dE)*Tkpsis(:,2);
    for chebkunt=2:kunt-1
        chebsum = chebsum + 2*((-1i)^chebkunt)*besselj(chebkunt,dt*dE)*Tkpsis(:,chebkunt+1);
    end
    
    psis(:,nrn) = chebsum;
    
    %% ploting
    
    plot(x,V,'color','g','linewidth',2);
    hold on
    plot(x,real(chebsum),'--','color','r')
    plot(x,imag(chebsum),'--','color','b')
    pdf = chebsum.*conj(chebsum);
    plot(x,pdf,'-','linewidth',1,'color','k','linewidth',2);
    axis([min(x), max(x), -1, 4]);
    %axis off
    hold off
    F = getframe(fig1);
    aviobj = addframe(aviobj,F);
    
    % record the autocorrelation...
    Pt(nrn) = trapz(x,chebsum.*conj(psi0));
    trapz(x,psis(:,nrn).*conj(psis(:,nrn)))

end
close(fig1);
aviobj = close(aviobj);


%% Form PE and get the spectrum...
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
plot(E,abs(Pe));

%% Practice points for the L1 reconstruction.

%seed the random generator with the example seeder
stream = RandStream('mrg32k3a');

num_samples = Nt/4;

%make our downsampling matrix
sample_points = randsample(stream,1:Nt, num_samples);

% FPC_AS A_operator class
A = A_operator( @(z) pifft(z,sample_points), @(z) pfft(z,sample_points,n) );

% tiny mu corresponds to heavy weight on the fidelity term
mu = 1e-10;

%downsample.
%u_samples are what we measure from the expensive machine
%we simulate measured data with the following line:
%u_samples = Po(sample_points)'; %apply the first operator in A_operator.
u_samples = transpose(Po(sample_points)*sqrt(length(Po)));

%fprintf('u-Po %f\n',norm(u - Po));

% Call Wotao's code.
%[Pe_new, Out] = FPC_AS(n, A, u_samples, mu);

[uHat_approx, Out] = FPC_AS(Nt, A, u_samples, mu);

close all;
figure()
hold on;
plot(abs(fft(usparseIFT)/sqrt(length(usparse)) ));
plot(abs(uHat_approx),'r--');

%% Neat picture
figure();
imshow(abs(psis(min(find(-L/9<x)):max(find(x<L/9)), :)))
colormap hot;
