clear all; close all; clc;

Ae = 4.0;

N = [2:12]';
error = zeros(size(N,1),1);
for i=1:size(N,1)
    [z,w] = zwgll(N(i));
    Bh = spdiags(w,0,N(i)+1,N(i)+1);
    F = ones(N(i)+1,N(i)+1);

    A = Bh*F*Bh';
    A = sum(sum(A));
    error(i) = abs(A-Ae);
end

semilogy(N,error/Ae); 
xlabel('N'); ylabel('Relative Error'); title('Square Domain');
axis square;