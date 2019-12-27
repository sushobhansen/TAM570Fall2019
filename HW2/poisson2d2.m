function U = poisson2d2(N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[z,w] = zwgll(N);
Bh = spdiags(w,0,N+1,N+1);
Dh = deriv_mat(z);
R = speye(N+1); R = R(2:end-1,:);

Ah = Dh'*Bh*Dh; Ah = 0.5*(Ah+Ah');
Ah = R*Ah*R';
Bhr = R*Bh*R';

[Sh,Lh] = eigs(Ah,Bhr,N-1);

Sh = bsxfun(@rdivide,Sh,vecnorm(Sh'*Bhr*Sh));
F = 1.0*exp(0*Bh);
Ih = speye(N-1);
L = kron(Ih,Lh)+kron(Lh,Ih);

U = (Sh'*((R*Bh)*F*(Bh'*R'))*Sh);
%U = L^(-1)*reshape(U,(N-1)^2,1);
%U = Sh*reshape(U,N-1,N-1)*Sh';
U = reshape(diag(L).^(-1),N-1,N-1).*U; %L is diagonal, this is cheap
U = Sh*U*Sh';
%U = padarray(U,[1,1],0,'both');
U = R'*U*R;

% [X,Y] = ndgrid(z,z);
% surf(X,Y,U);
% xlabel('x'); ylabel('y'); zlabel('u');
% title(sprintf('2D Poisson Eq for N = %d',N));
% axis square;
end

