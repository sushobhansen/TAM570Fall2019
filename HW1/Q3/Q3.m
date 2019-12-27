clear all; clc;

Nmax = 50;
etime = zeros(Nmax,1);
iters = 1;

for N=1:Nmax
    N
    M = ceil(1.5*N);
    A = randn(M+1,N+1);
    B = randn(M+1,N+1);
    C = randn(M+1,N+1);
    U = randn(N+1,N+1,N+1);
    
    for L=1:iters
        tstart = tic;
        V = kron3d1(A,B,C,U);
        etime(N) = etime(N) + toc(tstart);
    end
    
end

loglog(1:Nmax,etime)

etime = etime/iters;

