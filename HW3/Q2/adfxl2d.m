clear all; close all; clc; format compact; format longe;

Nx =41;

lx = 2*pi;
hx = lx/Nx; 
x = hx*[0:Nx-1]';
Mx = ceil(1/.5*(Nx+1)); if mod(Mx,2)==0; Mx=Mx+1; end;
hxm = lx/Mx; xm = hxm*[0:Mx-1]';
Df = (2*pi/lx)*dhatf(Nx);
Dff = Df^2;
Bhx = (0.5*lx/pi)*pi*speye(Nx); Bhx(1,1) = (0.5*lx/pi)*2*pi;
Bmx = (0.5*lx/pi)*pi*speye(Mx); Bmx(1,1) = (0.5*lx/pi)*2*pi;
P = sparse([eye(Nx);zeros(Mx-Nx,Nx)]);
Ahx = -Dff*Bhx;

Ny = 20;
ly = 2;
[y,wy] = zwgll(Ny);
Bhy = (ly/2)*spdiags(wy,0,Ny+1,Ny+1);
My = ceil(1.5*(Ny+1));
[ym,wym] = zwgll(My);
Jy = interp_mat(ym,y);
Bmy = (ly/2)*spdiags(wym,0,My+1,My+1);
Dhy = (2/ly)*deriv_mat(y);
Ry = eye(Ny+1); Ry = Ry(2:end-1,:);
Ahy = (Dhy'*Bhy*Dhy); Ahy = 0.5*(Ahy+Ahy');
Bhyr = Ry*Bhy*Ry';
[Sy,Ly] = eigs(Ry*Ahy*Ry',Bhyr,size(Ry,1));
Sy = bsxfun(@rdivide,Sy,vecnorm(Sy'*Bhyr*Sy));

[X,Y] = ndgrid(x,y);
[Xm,Ym] = ndgrid(xm,ym);

U0 = zeros(size(X));
U0(X>=1 & X<=2) = 1;

Cxm = 1-Ym.^2; Cxm=0*Cxm;
Cym = 0*Xm;  

c = max(max(max(abs(Cxm))),max(max(abs(Cym)))); 
dx = min(min(abs(diff(x))),min(abs(diff(y))));
nu = 0.001;
CFL = 0.6335/pi;
dt = CFL*dx/c;
tmax = 10;
nsteps = ceil(tmax/dt); 
nsteps = 1000;
dt = tmax/nsteps;
norm2 = zeros(nsteps,1);
norminf = norm2;

U = rfft(U0);
U = U*Ry';

Ihx = speye(Nx);
Ihy = speye(size(Ry,1));
for i=1:nsteps
    if i==1
        U0=U;
        [alpha,beta] = timecoeffs(i);
        D = alpha(1)*kron(Ihy,Bhx) + nu*dt*(kron(Ly,Bhx)+kron(Ihy,Ahx));
        D = reshape(diag(D),Nx,size(Ry,1));
        U = -alpha(2)*(Bhx*U0*Ry*Bhy'*Ry')-beta(1)*dt*P'*(Bmx*rfft(Cxm.*irfft(P*Df*U0*Ry*Jy'))*Bmy')*Jy*Ry';
        U = ((D.^(-1)).*(U*Sy))*Sy';
        U1=U;
    elseif i==2
        [alpha,beta] = timecoeffs(i);
        D = alpha(1)*kron(Ihy,Bhx) + nu*dt*(kron(Ly,Bhx)+kron(Ihy,Ahx));
        D = reshape(diag(D),Nx,size(Ry,1));
        U = -alpha(3)*(Bhx*U0*Ry*Bhy'*Ry')-alpha(2)*(Bhx*U1*Ry*Bhy'*Ry');
        U = U -beta(1)*dt*P'*(Bmx*rfft(Cxm.*irfft(P*Df*U1*Ry*Jy'))*Bmy')*Jy*Ry' - beta(2)*dt*P'*(Bmx*rfft(Cxm.*irfft(P*Df*U0*Ry*Jy'))*Bmy')*Jy*Ry';
        U = ((D.^(-1)).*(U*Sy))*Sy';
        U2=U;
    else
        [alpha,beta] = timecoeffs(3);
        D = alpha(1)*kron(Ihy,Bhx) + nu*dt*(kron(Ly,Bhx)+kron(Ihy,Ahx));
        D = reshape(diag(D),Nx,size(Ry,1));
        U = -alpha(4)*(Bhx*U0*Ry*Bhy'*Ry')-alpha(3)*(Bhx*U1*Ry*Bhy'*Ry')-alpha(2)*(Bhx*U2*Ry*Bhy'*Ry');
        U = U -beta(1)*dt*P'*(Bmx*rfft(Cxm.*irfft(P*Df*U2*Ry*Jy'))*Bmy')*Jy*Ry' - beta(2)*dt*P'*(Bmx*rfft(Cxm.*irfft(P*Df*U1*Ry*Jy'))*Bmy')*Jy*Ry' - beta(3)*dt*P'*(Bmx*rfft(Cxm.*irfft(P*Df*U0*Ry*Jy'))*Bmy')*Jy*Ry';
        U = ((D.^(-1)).*(U*Sy))*Sy';
        U0=U1; 
        U1=U2;
        U2=U;        
    end
    
    norm2(i) = sqrt(sumsqr(irfft(U*Ry))/((Ny+1)*Nx));
    norminf(i) = max(max(abs(irfft(U*Ry))));
    
end
 

Ub = irfft(U*Ry);
mesh(X,Y,Ub); %colorbar;
title(sprintf('Time = %f, step=%d',tmax,nsteps));
figure;
semilogy(1:nsteps,norm2,1:nsteps,norminf); xlabel('Time step'); ylabel('norm');
legend('2-Norm','Inf Norm');