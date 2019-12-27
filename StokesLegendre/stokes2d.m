clear all; close all; clc; format compact; format longe;
N = 8; k=2; Np = N-k;
[z,w] = zwgll(N);
[zp,wp] = zwgll(Np);
Bh = spdiags(w,0,N+1,N+1);
B = kron(Bh,Bh);
Jh = interp_mat(z,zp);
J = kron(Jh,Jh);
Dh = deriv_mat(z);
Rx = speye(N+1); Rx = Rx(2:end-1,:);
Ry = speye(N+1); Ry = Ry(2:end-1,:);
R1 = kron(Ry,Rx);
Ahx = Dh'*Bh*Dh;
Ahy = Dh'*Bh*Dh;
A = kron(Bh,Ahx)+kron(Ahy,Bh);
A1 = R1*A*R1'; %Matrix for u1

Rx = speye(N+1); Rx = Rx(2:end,:);
Ry = speye(N+1); Ry = Ry(2:end-1,:);
R2 = kron(Ry,Rx);
Ahx = Dh'*Bh*Dh;
Ahy = Dh'*Bh*Dh;
A = kron(Bh,Ahx)+kron(Ahy,Bh);
A2 = R2*A*R2'; %matrix for u2

Ih = speye(N+1);
AA = [A1 zeros(size(A1,1),size(A2,2)); zeros(size(A2,1),size(A1,2)) A2];
Dx = kron(Ih,Dh)*R1';
Dy = kron(Dh,Ih)*R2';
DD = [J'*B*Dx J'*B*Dy];

S = DD*AA^(-1)*DD';
[V,e] = eig(S); e = sort(abs(diag(e))); e = e(1:20);
semilogy(1:20,e,'ro'); xlabel('Eigenvalue number'); ylabel('|Eigenvalue|');
axis square; grid minor; xlim([0 20]);
title(sprintf('N = %d, N_p = %d, Symmetry on x=1',N,Np));
