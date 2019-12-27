clear all; close all; clc;

%% Mesh
Nmax = 15;
niters = 0*[2:Nmax]';
ntimes = niters;
niters_precon = niters;
ntimes_precon = niters;
for N=2:Nmax
    N
    [z,w] = zwgll(N);
    [X,Y,Z] = mesh3d(N);
    lx = (sqrt(2)+1)/2;
    ly = 1/(2*sqrt(2));
    lz = 0.1;

    %% Evaluate derivatives 
    Dh = deriv_mat(z);
    Bh = spdiags(w,0,N+1,N+1);
    Ih = speye(N+1);

    Xr = kron3d2(Ih,Ih,Dh,X); 
    Xs = kron3d2(Ih,Dh,Ih,X);
    Xt = kron3d2(Dh,Ih,Ih,X);

    Yr = kron3d2(Ih,Ih,Dh,Y); 
    Ys = kron3d2(Ih,Dh,Ih,Y);
    Yt = kron3d2(Dh,Ih,Ih,Y);

    Zr = kron3d2(Ih,Ih,Dh,Z); 
    Zs = kron3d2(Ih,Dh,Ih,Z);
    Zt = kron3d2(Dh,Ih,Ih,Z);

    J = Xr.*(Ys.*Zt - Yt.*Zs) - Xs.*(Yr.*Zt - Yt.*Zr) + Xt.*(Yr.*Zs - Ys.*Zr);
    J = 1./J;

    rx = (Ys.*Zt-Zs.*Yt).*J; ry = (Xt.*Zs - Xs.*Zt).*J;  rz = (Xs.*Yt-Ys.*Xt).*J;
    sx = (Yt.*Zr-Yr.*Zt).*J; sy = (Xr.*Zt-Xt.*Zr).*J; sz = (Yr.*Xt-Yt.*Xr).*J;
    tx = (Yr.*Zs-Ys.*Zr).*J; ty = (Xs.*Zr-Xr.*Zs).*J; tz = (Xr.*Ys-Xs.*Yr).*J;

    J = 1./J;

    %% Define restriction matrices
    %Rx = Ih;
    Rx = Ih(2:end-1,:); %Verification
    %Ry = Ih(2:end,:);
    Ry = Ih(2:end-1,:); %Verification
    %Rz = Ih;
    Rz = Ih(2:end-1,:); %Verification
    %F = exp(10*Y.*Z).*sin(10*X); 
    F = 1+0*X; %For verification
    B = kron(Bh,kron(Bh,Bh));

    %% Define RHS
    b = kron3d2(Rz*Bh,Ry*Bh,Rx*Bh,J.*F);

    %%Solve and save values
    preconflag = 0;
    tstart = tic;
    [U,flag,relres,iter,resvec] = solve(J,B,Dh,Rx,Ry,Rz,rx,ry,rz,sx,sy,sz,tx,ty,tz,b,lx,ly,lz,N,preconflag);
    ntimes(N-1) = toc(tstart);
    niters(N-1) = iter;
    lastiter = iter;
    
    preconflag = 1;
    tstart = tic;
    [U,flag,relres,iter,resvec_precon] = solve(J,B,Dh,Rx,Ry,Rz,rx,ry,rz,sx,sy,sz,tx,ty,tz,b,lx,ly,lz,N,preconflag);
    ntimes_precon(N-1) = toc(tstart);
    niters_precon(N-1) = iter;
    lastiter_precon = iter;
    
    pause(0.01);
end

semilogy(2:Nmax,niters,'r.-',2:Nmax,niters_precon,'b-o',2:Nmax,10*(2:Nmax).^1.5,'k');
xlabel('N'); ylabel('Number of iterations'); title('Iterations to Convergence');
legend('CG','PCG','O(N^{1.5})','Location','northwest');
axis square;
figure;
semilogy(2:Nmax,ntimes,'r.-',2:Nmax,ntimes_precon,'b-o',2:Nmax,10*(2:Nmax).^1.5,'k');
xlabel('N'); ylabel('Solve time (s)'); title('Convergence Time');
legend('CG','PCG','O(N^{1.5})','Location','northwest');
axis square;
figure;
semilogy(0:lastiter,resvec,'r-',0:lastiter_precon,resvec_precon,'b',0:lastiter,(0:lastiter).^-2,'k',0:lastiter,(0:lastiter).^-3,'g');
xlabel('Iteration'); ylabel('residual'); title('Residuals at N=30');
legend('CG','PCG','O(N^2)','O(N^4)');
axis square;
figure;
%% Plot
U = reshape(U,size(Rz,1),size(Ry,1),size(Rx,1));
U = kron3d2(Rz',Ry',Rx',U);
surf(X(:,:,end-1),Y(:,:,end-1),U(:,:,end-1)); 
xlabel('x'); ylabel('y'); zlabel('u'); 
title('f \equiv 0, deformed geometry, u(\partial \Omega)=0');
axis square;