function pot = V(r)

De = 0.01688;
re = 3.28892;
b = 1.47612;


pot = De*(1-exp(-b*(r-re)))^2 - De;