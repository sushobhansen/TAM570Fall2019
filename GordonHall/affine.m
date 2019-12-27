clear all; close all; clc;
Ae = 0.25;
N = [2:12]';
error = zeros(size(N,1),1);
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

    Xr = Dh*X; Xs = X*Dh';
    Yr = Dh*Y; Ys = Y*Dh';
    J = (Xr.*Ys - Yr.*Xs);
    %mesh(X,Y,0*X)

    F = ones(N(i)+1,N(i)+1);
    A = sum(sum(Bh*(J.*F)*Bh'));
    error(i) = abs(A-Ae);
end

semilogy(N,error/Ae);
xlabel('N'); ylabel('Relative Error'); title('Affine Transformation');
axis square;