clear all; clc;

N = 9; M = 50;
[z,w] = zwgll(N);
[zf,wf] = zwuni(M);

Jh = interp_mat(zf,z);

n = N+1;
U = zeros(n,n); U(1,1) = 1;
Uf = Jh*U*Jh';

[X,Y] = ndgrid(z,z);
Ux = Jh*U;
Xx = Jh*X; Yx = Jh*Y;

Xf = Jh*X*Jh'; 
Yf = Jh*Y*Jh';
Uxf = Jh*U*Jh';
surf(Xf,Yf,Uf,'EdgeColor','none'); hold on;
xlabel('r'); ylabel('s'); zlabel('U');
plot3(Xf,Yf,Uxf); 
plot3(Xf',Yf',Uxf');

Xy = X*Jh';
Yy = Y*Jh';
Uy = U*Jh';
%plot3(Xy,Yy,Uy)

%% Distorted geometry
Xc = [0 1;0 2]'; v
Yc = [0 0;1 2]';
[z1,w1] = zwgll(1);
Jz = interp_mat(z,z1);
X = Jz*Xc*Jz'; Y = Jz*Yc*Jz';
figure;
%surf(X,Y,U)

Xf = Jh*X*Jh'; Yf = Jh*Y*Jh';
surf(Xf,Yf,Uf,'EdgeColor','none'); hold on;
plot3(Xf,Yf,Uxf);
plot3(Xf',Yf',Uxf');
xlabel('x'); ylabel('y'); zlabel('U');