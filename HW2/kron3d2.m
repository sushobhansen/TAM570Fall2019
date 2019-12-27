function V = kron3d2(A,B,C,U)
%Smart implementation of Kronicker product
    
    N = size(A,2); N=N-1;
    M = size(A,1); M=M-1;
    V = zeros(M+1,M+1,M+1);
    F = zeros(M+1,M+1,N+1);
    
    U = permute(U,[1,3,2]);
    
    for r = 1:N+1
        F(:,:,r) = C*U(:,:,r)*A';
    end
    
    F  = permute(F,[3,1,2]);
    for j = 1:M+1
        V(:,:,j) = V(:,:,j)+B*F(:,:,j);
    end
    
    V = permute(V,[2,1,3]);
end

