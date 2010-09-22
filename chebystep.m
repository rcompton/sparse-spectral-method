function chebsum = chebystep(M,dx,L,V,dt,psi0,plot)
%
% Apply a Chebyshev polynomial approximation
% to propagate the exponential by dt
%
n = max(size(V));

dE = (pi^2)/(2*M*dx^2) + max(V) - min(V);
maxE = (pi^2)/(2*M*dx^2) + max(V);
minE = min(V);

%Hnorm = 2*(H - minE)/(maxE - minE) - eye(n);


%figure maxk for a given timestep
maxk = 50;
nexttenjays = 1:.3:1.5;
while max( abs(besselj(maxk,dE*dt*nexttenjays)) ) > 1e-6
    maxk = maxk+1; 
end
maxk;
dE*dt;

k = -n/2:(n/2-1);
k = k.*pi./L;
k = fftshift(k);

Hnormspec = inline( '(2/dE)*((-1/(2*M))*ifft(-(k.^2).*fft(Phi)) + V.*Phi) - (1+2*minE/dE).*Phi','M', 'V','dE','minE','k','dx','Phi');


%make maxk Tkpsis
Tkpsis = zeros(maxk,n);
Tkpsis(1,:) = psi0;

%Hu = (-1/(2*M))*([0 0 psi0] - 2*[0 psi0 0] + [psi0 0 0])/(dx^2); % Hnorm*psi0;
%Tkpsis(2,:) = Hu(2:end-1) + V.*psi0;

Tkpsis(2,:) = Hnormspec(M,V,dE,minE,k,dx,psi0);

for k=3:maxk
    %Hu = (-1/(2*M))*([0 0 Tkpsis(k-1,:)] - 2*[0 Tkpsis(k-1,:) 0] + [Tkpsis(k-1,:) 0 0])/(dx^2); % Hnorm*psi0;
    %Hu = Hu(2:end-1) + V.*Tkpsis(k-1,:);
    %Hu = (2/dE)*Hu - ((2*minE/dE) + 1)*Tkpsis(k-1,:);
    Hu = Hnormspec(M,V,dE,minE,k,dx,Tkpsis(k-1,:) );
    Tkpsis(k,:) = 2*Hu - Tkpsis(k-2,:);
end

%make the sum
chebsum = besselj(0,dt*dE)*Tkpsis(1,:);
chebsum = chebsum + -2i*besselj(1,dt*dE)*Tkpsis(2,:);
for k=2:maxk-1
    chebsum = chebsum + 2*((-1i)^k)*besselj(k,dt*dE)*Tkpsis(k+1,:);
    %bsls(k) = abs(besselj(k,dt*dE));
end

% if plot
%     figure();
%     plot(bsls);
% end
