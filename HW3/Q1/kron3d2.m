function V = kron3d2(A,B,C,U)
%Smart implementation of Kronicker product
    
    Na = size(A,2); Na=Na-1;
    Ma = size(A,1); Ma=Ma-1;
    
    Nb = size(B,2); Nb=Nb-1;
    Mb = size(B,1); Mb=Mb-1;
    
    Nc = size(C,2); Nc=Nc-1;
    Mc = size(C,1); Mc=Mc-1;
    
    Nu = size(U,2); Nu=Nu-1;
    Mu = size(U,1); Mu=Mu-1;
    
    V = zeros(Mb+1,Ma+1,Mc+1);
    F = zeros(Mc+1,Ma+1,Nu+1);
    
    U = permute(U,[1,3,2]);
    
    for r = 1:Nu+1
        F(:,:,r) = C*U(:,:,r)*A';
    end
    
    F  = permute(F,[3,1,2]);
    for j = 1:Mc+1
        V(:,:,j) = V(:,:,j)+B*F(:,:,j);
    end
    
    V = permute(V,[2,1,3]);
end

