clear all; close all; clc;
N = [60]';
errors = zeros(size(N,1),1);

for j=1:size(N,1)
    N(j)
    M = ceil(3*(N(j)+1)/2);

    [z,w] = zwgll(N(j));
    Bh = spdiags(w,0,N(j)+1,N(j)+1);
    [zm,wm] = zwgll(M);
    Bmh = spdiags(wm,0,M+1,M+1);
    R = speye(N(j)+1); R = R(2:end-1,:);
    Jh = interp_mat(zm,z);
    Dh = deriv_mat(z);
    Dtilde = Jh*Dh;
    [Xm,Ym] = ndgrid(zm,zm);
    Cxm = -Ym; Cym = Xm;

    Bm = kron(Bmh,Bmh); %Sparse matrix, cheap

    h = min(diff(z));
    CFL = 1.0;
    c = 1.0;
    dt = CFL*h/c;
    T = 2*pi;
    steps = ceil(T/dt);

    [X,Y] = ndgrid(z,z);
    x0 = 0.5; y0 = 0.0;
    r2 = (X-x0).^2+(Y-y0).^2;
    U0 = exp(-r2/0.016); %IC
    U = U0;
    for i=1:steps
        %surf(X,Y,U);
        %title(spriN(j)tf('Time = %f',dt*(i-1)));
        %pause(0.01);

        %RK4    
        k1 = dt*advect2drhs(Bh,Bm,Cxm,Cym,Dtilde,Jh,R,U);
        k2 = dt*advect2drhs(Bh,Bm,Cxm,Cym,Dtilde,Jh,R,U+k1/2);
        k3 = dt*advect2drhs(Bh,Bm,Cxm,Cym,Dtilde,Jh,R,U+k2/2);
        k4 = dt*advect2drhs(Bh,Bm,Cxm,Cym,Dtilde,Jh,R,U+k3);
        U = U + (k1+2*k2+2*k3+k4)/6.0;
    end
    
    %errors(j) = max(max(abs(U-U0)));
    errors(j) = abs(max(max(abs(U)))-max(max(abs(U0))));
end

% semilogy(N,errors,N,1000*N.^(-4)); xlabel('N'); ylabel('$||u^n-u^0||_\infty$');
% semilogy(N,errors,N,1000*N.^(-4)); xlabel('N'); ylabel('||u^n-u^0||_\infty'); title(sprintf('CFL = %f',CFL));
% axis square;
% legend('Error','O(N^4)');

surf(X,Y,U);
xlabel('x'); ylabel('y'); zlabel('u');
title(sprintf('CFL = %f, N = %d',CFL,N(1)));
axis square;
