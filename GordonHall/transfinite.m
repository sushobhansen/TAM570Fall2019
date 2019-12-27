clear all; close all; clc;

N = [2:12]';
error = zeros(size(N,1),1);
Ae = 0.25*(pi-1.0);
for i=1:size(N,1)
    [z,w] = zwgll(N(i));
    [z1,w1] = zwgll(1);
    J1 = interp_mat(z,z1);
    Dh = deriv_mat(z);
    Bh = spdiags(w,0,N(i)+1,N(i)+1);
    X1 = [-1/2 1/2;-1/sqrt(2) 1/sqrt(2)]';
    Y1 = [1/2 1/2;1/sqrt(2) 1/sqrt(2)]';
    X = J1*X1*J1';
    Y = J1*Y1*J1';

    R = 1.0;
    x0 = 0.0; y0 = 0.0;
    delta = y0 + sqrt(R^2 - (X(:,end)-x0).^2) - Y(:,end);
    Y = Y + delta*J1(:,end)';
    %mesh(X,Y,0*Y);

    Xr = Dh*X; Xs = X*Dh';
    Yr = Dh*Y; Ys = Y*Dh';
    J = (Xr.*Ys - Yr.*Xs);
    
    F = ones(N(i)+1,N(i)+1);
    A = sum(sum(Bh*(J.*F)*Bh'));
    error(i) = abs(A-Ae);
end

semilogy(N,error/Ae); 
xlabel('N'); ylabel('Relative Error'); title('Transfinite Interpolation');
axis square;