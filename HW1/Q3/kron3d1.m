function V = kron3d1(A,B,C,U)
%Naive implementation of Kronicker product


    N = size(A,2); N=N-1;
    M = size(A,1); M=M-1;
    U = reshape(U,[(N+1)^3,1]);
    V = kron(kron(A,B),C)*U;
    V = reshape(V,[M+1,M+1,M+1]);

end

