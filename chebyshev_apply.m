function [Udt] = chebyshev_apply(dx, M, V, k, dt, Phi)
    %determine number of terms in Chebyshev expansion
    R = 0.5*dt*( (pi^2)/(2*M*dx^2) + max(V) - min(V) );
    G = min(V)*dt;
    
    %Like in Kosloff's orig paper
    alpha = 1.5;
    nterms = ceil(alpha*R)
    nterms = 100;
    
    R
    
    X = inline( '-1i*(dt/R)*(ifft(-(k.^2).*fft(Phi)) + V.*Phi )','V','k','R','dt','Phi'); 
    
    Phikm2 = Phi;
    Phikm1 = X(V,k,R,dt,Phi);
    
    Udt = 0 + 1*besselj(1,R)*Phikm1;
    
    %recursion
    count = 2;
    while besselj(count, R) > 1e-7
        Phik = 2*X(V,k,R,dt,Phikm1) + Phikm2;
        Phikm2 = Phikm1;
        Phikm1 = Phik;
        
        Udt = Udt + 2*besselj(count, R).*Phik;
        
        count = count+1;
        
        if mod(count, 100) == 1
            count
            besselk(count,R)
        end
        
    end