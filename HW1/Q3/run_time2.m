for L = 1:5
ftime = zeros(1,128);    
for N = 1:128
    N
tstart = tic;
M = ceil(3*N/2);
u = randn(N+1,N+1,N+1);
u = permute(u,[1,3,2]);
A = randn(M+1,N+1);
B = randn(M+1,N+1);
C = randn(M+1,N+1);
V = zeros(M+1,M+1,M+1);
F = zeros(M+1,M+1,N+1);
%W = zeros((M+1)^3,1);
for r = 1:N+1
    F(:,:,r) = C*u(:,:,r)*A';
end
F  = permute(F,[3,1,2]);
for j = 1:M+1
    V(:,:,j) = V(:,:,j)+B*F(:,:,j);
end
V = permute(V,[2,1,3]);
etime(N) = toc(tstart);
end
ftime = etime + ftime;
end
ftime = ftime/L;
plot(1:128,ftime)