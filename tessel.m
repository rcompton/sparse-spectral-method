%
% Look at Bessel functions
% study of the decay as k exceeds R
% when k>R J_k(R) drops to 0 almost right away
% when k<R J_k(R) is not worth studying
%
clear all;close all;clc;

ks = 1:300;
Rs = 1:.5:500;

m = length(ks);
n = length(Rs);

bessels = zeros(m,n);
for i=1:m
    for j=1:n
        bessels(i,j) = besselj(ks(i),Rs(j));
    end
end

%figure();
%plot3(ks,Rs,bessels);

figure();
imshow(bessels); colormap hot;