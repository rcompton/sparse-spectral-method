N = 1000;
x = linspace(-1,1,N);


t = 20;

trapz(x,exp(-1i*x*t).*chebpoly(8,x)./sqrt(1-x.^2))

besselj(8,t)