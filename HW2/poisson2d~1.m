function U = poisson2d(N,M)
%Solves the 2D Poisson Eq with homogeneous Dirichlet conditions
[z,w] = zwgll(N);
[zm,wm] = zwgll(M);
Bh = spdiags(w,0,N+1,N+1);
Bm = spdiags(wm,0,M+1,M+1);
Dh = deriv_mat(z);
R = speye(N+1); R = R(2:end-1,:);
Rm = speye(M+1); Rm = Rm(2:end-1,:);
Jh = interp_mat(zm,z);

Ah = Dh'*Jh'*Bm*Jh*Dh; Ah = 0.5*(Ah+Ah');
Ah = R*Ah*R';
Bh = R*Bh*R';
Bmj = R*Jh'*Bm*Jh*R';
Bm = Rm*Bm*Rm';
Jh = Rm*Jh*R';

[Sh,Lh] = eigs(Ah,Bmj,N-1);

Sh = bsxfun(@rdivide,Sh,vecnorm(Sh'*Bmj*Sh));
F = 1.0*exp(0*Bm);
Ih = spdiags(diag(Sh'*Bmj*Sh),0,N-1,N-1);
L = kron(Ih,Lh)+kron(Lh,Ih);

U = (Sh'*(Jh'*Bm*F*Bm'*Jh)*Sh);
U = L^(-1)*reshape(U,(N-1)^2,1); %L is diagonal, this is cheap
U = Sh*reshape(U,N-1,N-1)*Sh';
U = padarray(U,[1,1],0,'both');

%[X,Y] = ndgrid(z,z);
%surf(X,Y,U)

end

