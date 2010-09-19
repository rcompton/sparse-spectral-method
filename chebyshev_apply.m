function [Udt] = chebyshev_apply(dx, M, V, k, dt, Phi)
    R = 0.5*dt*( (pi^2)/(2*M*dx^2) + max(V) - min(V) );
    G = min(V)*dt;
            
    
    %Like in Kosloff's orig paper
    %alpha = 1.5;
    %nterms = ceil(alpha*R)
    %nterms = 100;
    %R
    
    %X = inline( '-1i*(dt/R)*(ifft(-(k.^2).*fft(Phi)) + V.*Phi )','V','k','R','dt','Phi'); 
    
    %HPhi = (-1/(2*M))*([0 0 Phi] - 2*[0 Phi 0] + [Phi 0 0])/(dx^2);
    %HPhi = HPhi(2:end-1) + V.*Phi;
    
    %X = inline( '-1i*(dt/R)*([0 0 Phi] - 2*[0 Phi 0] + [Phi 0 0])/(dx^2) 
    X = inline('-1i*(dt/R)*( (-1/(2*M))*del2(Phi,dx) + V.*Phi) ', 'V','R','dx','dt','M','Phi');
    
    Phikm2 = Phi;
    %Phikm1 = X(V,k,R,dt,Phi);
    Phikm1 = X(V,R,dx,dt,M,Phi);
    
    Udt = 0 + 1*besselj(1,R)*Phikm1;
    
    
    %recursion, truncate when coefficient shrinks enough.
    count = 2;
    while besselj(count, R) > 1e-9
        %Phik = 2*X(V,k,R,dt,Phikm1) + Phikm2;
        Phik = 2*X(V,R,dx,dt,M,Phikm1) + Phikm2;
        Phikm2 = Phikm1;
        Phikm1 = Phik;
        
        
        
        Udt = Udt + 2*besselj(count, R).*Phik;
        
        count = count+1
        besselj(count,R)
        
        if mod(count, 10) == 1
            count
            besselj(count,R)
        end
        
    end