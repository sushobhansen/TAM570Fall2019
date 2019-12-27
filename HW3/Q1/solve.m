function [U,flag,relres,iter,resvec] = solve(J,B,Dh,Rx,Ry,Rz,rx,ry,rz,sx,sy,sz,tx,ty,tz,b,lx,ly,lz,N,preconflag)

b = reshape(b,size(Rz,1)*size(Ry,1)*size(Rx,1),1);
if(preconflag)
    [U,flag,relres,iter,resvec] = pcg(@multiplier,b,1e-10,10000,@preconditioner);
else
    [U,flag,relres,iter,resvec] = pcg(@multiplier,b,1e-10,10000);
end

% n = 3*3*2;
% w = zeros(n,n);
% for j=1:n
%    u = zeros(n,1); u(j)=1;
%    w(:,j) = multiplier(u);
% end
% U = w;

    function conditionedvec = preconditioner(U)
        [z,w] = zwgll(N);
        Bh = spdiags(w,0,N+1,N+1);
        Dh = deriv_mat(z);
        U = reshape(U,size(Rx,1),size(Ry,1),size(Rz,1));

        Ah = Dh'*Bh*Dh; Ah = 0.5*(Ah+Ah');
        Ahx = (2/lx)*Ah; Ahy = (2/ly)*Ah; Ahz = (2/lz)*Ah;
        Bhx = (lx/2)*Bh; Bhy = (ly/2)*Bh; Bhz = (2/lz)*Bh;
        Ax = Rx*Ahx*Rx'; Ay = Ry*Ahy*Ry'; Az = Rz*Ahz*Rz';
        Bx = Rx*Bhx*Rx'; By = Ry*Bhy*Ry'; Bz = Rz*Bhz*Rz';

        [Sx,Lx] = eigs(Ax,Bx,size(Ax,1));
        Sx = bsxfun(@rdivide,Sx,vecnorm(Sx'*Bx*Sx));

        [Sy,Ly] = eigs(Ay,By,size(Ay,1));
        Sy = bsxfun(@rdivide,Sy,vecnorm(Sy'*By*Sy));

        [Sz,Lz] = eigs(Az,Bz,size(Az,1));
        Sz = bsxfun(@rdivide,Sz,vecnorm(Sz'*Bz*Sz));

        Ix = speye(size(Ax,1)); Iy = speye(size(Ay,1)); Iz = speye(size(Az,1));
        L = kron(Iz,kron(Iy,Lx)) + kron(Iz,kron(Ly,Ix)) + kron(Lz,kron(Iy,Ix));

        conditionedvec = kron3d2(Sz',Sy',Sx',U);
        conditionedvec = reshape(full(diag(L).^(-1)),size(Rx,1),size(Ry,1),size(Rz,1)).*conditionedvec;
        conditionedvec = kron3d2(Sz,Sy,Sx,conditionedvec);

        conditionedvec = reshape(conditionedvec,size(Rx,1)*size(Ry,1)*size(Rz,1),1);
    end

    function integral = multiplier(U)
        N = size(Dh,1); N=N-1;
        Bdiag = reshape(full(diag(B)),N+1,N+1,N+1);
        U = reshape(U,size(Rx,1),size(Ry,1),size(Rz,1));
        Ix1 = rx.*kron3d2(Rz',Ry',Dh*Rx',U) + sx.*kron3d2(Rz',Dh*Ry',Rx',U) + tx.*kron3d2(Dh*Rz',Ry',Rx',U);
        Ix1 = Bdiag.*J.*Ix1;

        Iy1 = ry.*kron3d2(Rz',Ry',Dh*Rx',U) + sy.*kron3d2(Rz',Dh*Ry',Rx',U) + ty.*kron3d2(Dh*Rz',Ry',Rx',U);
        Iy1 = Bdiag.*J.*Iy1;
        
        Iz1 = rz.*kron3d2(Rz',Ry',Dh*Rx',U) + sz.*kron3d2(Rz',Dh*Ry',Rx',U) + tz.*kron3d2(Dh*Rz',Ry',Rx',U);
        Iz1 = Bdiag.*J.*Iz1;
        
        Ix = rx.*Ix1 + ry.*Iy1 + rz.*Iz1;
        Ix = kron3d2(Rz,Ry,Rx*Dh',Ix);
        
        Iy = sx.*Ix1 + sy.*Iy1 + sz.*Iz1;
        Iy = kron3d2(Rz,Ry*Dh',Rz,Iy);
        
        Iz = tx.*Iz1 + ty.*Iz1 + tz.*Iz1;
        Iz = kron3d2(Rz*Dh',Ry,Rz,Iz);
        
        integral = Ix+Iy+Iz;
        integral = reshape(integral,size(Rx,1)*size(Ry,1)*size(Rz,1),1);
    end

end

