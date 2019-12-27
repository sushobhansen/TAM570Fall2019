clear all; clc; close all;
N = 21;
xn = 2*pi*[0:N-1]'/N;
M = 10*N; 
xm = 2*pi*[0:M-1]'/M;

un = sign(xn-pi);
um = interp_f(M,un);
plot(xm,um,'r-',xn,un,'bo','linewidth',1.4);
hold on;
plot(xm,0*xm,'k','linewidth',2.0);

N = [20:2:200]';
maxum = zeros(size(N,1),1);
for i=1:size(N,1)
    xn = 2*pi*[0:N(i)-1]'/N(i);
    M = 10*N(i); 
    xm = 2*pi*[0:M-1]'/M;

    un = sign(xn-pi);
    um = interp_f(M,un);
    maxum(i) = max(abs(um));
end
figure;
plot(N,maxum);
