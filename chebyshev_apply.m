function [Udt] = chebyshev_apply(dx,M,V.dt,Phi)
    %determine number of terms in Chebyshev expansion
    R = 0.5*dt*( (pi^2)/(2*M*dx^2) + max(V) - min(V) );
    G = min(V)*dt;
    
    %Like in Kosloff's orig paper
    alpha = 1.5;
    nterms = alpha*R;
    
    X = inline( '-1i*(dt/R)*(ifft(-(k.^2).*fft(Phi)) + V.*Phi )','V','k','Phi'); 
    
    Phikm2 = Phi;
    Phikm1 = X(V,k,Phi);
    
    Udt = 0 + 1*besselj(1,R)*Phikm1;
    
    %recursion
    for i=3:nterms
        Phik = 2*X(Phikm1) + Phikm2;
        Phikm2 = Phikm1;
        Phikm1 = Phik;
        
        Udt = Udt + 2*besselj(k,R)*Phik;
        
    end