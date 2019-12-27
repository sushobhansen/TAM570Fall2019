clear all; close all; clc;
format compact; format longe;

N = 50;
U = poisson2d2(N);
Umaxe = max(max(abs(U)));

N = [2:2:40]';
error = zeros(size(N,1),1);
for i=1:size(N,1)
   U = poisson2d2(N(i)); 
   error(i) = abs(max(max(abs(U)))-Umaxe);
end

subplot(1,2,1)
loglog(N,error,N,N.^(-2),N,N.^(-4));
xlabel('N'); ylabel('error');
title('2D Poisson Convergence');
legend('Spectral','O(N^2)','O(N^4)');
axis square;

subplot(1,2,2)
semilogy(N,error,N,N.^(-2),N,N.^(-4));
xlabel('N'); ylabel('error');
title('2D Poisson Convergence');
legend('Spectral','O(N^2)','O(N^4)');
axis square;