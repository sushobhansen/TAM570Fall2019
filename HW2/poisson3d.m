function U = poisson3d(N)
%Solved Poisson 3D with D-N homogeneous BCs

[z,w] = zwgll(N);
Bh = spdiags(w,0,N+1,N+1);
Dh = deriv_mat(z);
R = speye(N+1); R = R(2:end,:);


Ah = Dh'*Bh*Dh; Ah = 0.5*(Ah+Ah');
Ah = R*Ah*R'; 
F = 1.0+0*Bh;
Bhr = R*Bh*R'; 

[Sh,Lh] = eigs(Ah,Bhr,size(Ah,1));
Sh = bsxfun(@rdivide,Sh,vecnorm(Sh'*Bhr*Sh));

Ih = speye(size(Ah,1));
L = kron(Ih,kron(Ih,Lh)) + kron(Ih,kron(Lh,Ih)) + kron(Lh,kron(Ih,Ih));

U = kron3d2(R*Bh,R*Bh,R*Bh,F);
U = kron3d2(Sh',Sh',Sh',F);
U = reshape(full(diag(L).^(-1)),N,N,N).*U;
U = kron3d2(Sh,Sh,Sh,reshape(U,N,N,N));

U = kron3d2(R',R',R',U);

%[X,Y] = ndgrid(z,z);
%surf(X,Y,Ub(:,:,2))

end

