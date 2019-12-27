function [X,Y,Z] = mesh3d(N)

[z,w] = zwgll(N);
[z1,w1] = zwgll(1);
J1 = interp_mat(z,z1);
Dh = deriv_mat(z);

X1 = [-1/2 1/2;-1/sqrt(2) 1/sqrt(2)]'; X1(:,:,2) = X1;
Y1 = [1/2 1/2;1/sqrt(2) 1/sqrt(2)]'; Y1(:,:,2) = Y1;
Z1 = [0 0;0 0]'; Z1(:,:,2) = ones(2,2)*0.1;

X = kron3d2(J1,J1,J1,X1);
Y = kron3d2(J1,J1,J1,Y1);
Z = kron3d2(J1,J1,J1,Z1);

R = 1.0;
x0 = 0.0; y0 = 0.0;
delta = y0 + sqrt(R^2 - (X(:,end,:)-x0).^2) - Y(:,end,:);

for i=1:N+1
Y(:,:,i) = Y(:,:,i) + delta(:,:,i)*J1(:,end)'; 
end

end

