clear all; close all; clc;
format longe; format compact;

ax=1; 
bx=6.0; 
ay=-1;
by=1.0;
lx = bx-ax; ly = by-ay;
N=20;
%Nx = 20; Ny = 20;
Nx=N; Ny=N;
Mx = ceil(1.5*(Nx+1)); My = ceil(1.5*(Ny+1));
Npx = Nx-2; Npy = Ny-2;
nu = 1e-4; g = 9.81; rho_0 = 1000; Pr = 1;
ntimes = 5;
comp_time = 0;

for n=1:ntimes
    tic;
%% U-mesh (P N) and P-mesh (P N-2)
[x,wx] = zwgll(Nx);
[y,wy] = zwgll(Ny);
[xm,wxm] = zwgll(Mx);
[ym,wym] = zwgll(My);
[xp,wpx] = zwgll(Npx);
[yp,wpy] = zwgll(Npy);

Jux = interp_mat(xm,x);
Juy = interp_mat(ym,y);
Jpx = interp_mat(x,xp);
Jpy = interp_mat(y,yp);

Bhx = (lx/2)*spdiags(wx,0,Nx+1,Nx+1);
Bhy = (ly/2)*spdiags(wy,0,Ny+1,Ny+1);
Bmx = (lx/2)*spdiags(wxm,0,Mx+1,Mx+1);
Bmy = (ly/2)*spdiags(wym,0,My+1,My+1);
Dhx = (2/lx)*deriv_mat(x);
Dhy = (2/ly)*deriv_mat(y);
Ahx = Dhx'*Bhx*Dhx; Ahx = 0.5*(Ahx+Ahx');
Ahy = Dhy'*Bhy*Dhy; Ahy = 0.5*(Ahy+Ahy');

x=ax+(x+1)*lx/2;
y=ay+(y+1)*ly/2;
xm=ax+(xm+1)*lx/2;
ym=ay+(ym+1)*ly/2;
xp=ax+(xp+1)*lx/2;
yp=ay+(yp+1)*ly/2;
[Xplot,Yplot]=meshgrid(x,y);
[X,Y] = ndgrid(x,y);

%% Boundaries
U1b = 0*X;
boundary = 0*X;
%Left  & Right
U1b(1,:) = boundary(1,:); U1b(end,:) = boundary(end,:);
%U1b(1,:) = 1-y.*y; U1b(end,:) = 1-y.*y;
%Down  &  Up
U1b(:,1) = boundary(:,1); U1b(:,end) = boundary(:,end);

U2b = 0*Y;

U2b(1,:) =boundary(1,:); U2b(end,:) = boundary(end,:);
%Down  &  Up
U2b(:,1) = boundary(:,1); U2b(:,end) = boundary(:,end);

Ihx = speye(Nx+1);
Ihy = speye(Ny+1);
Rx1 = Ihx(2:end-1,:); %D-D
Ry1 = Ihy(2:end-1,:); %D-D
%Ry1 = [speye(Ny) zeros(Ny,1)]; Ry1(1,1)=1; %Periodic

Rx2 = Ihx(2:end-1,:); %D-D
Ry2 = Ihy(2:end-1,:); %D-D
%Ry2 = [speye(Ny) zeros(Ny,1)]; Ry2(1,1)=1; %Periodic

rhobb = exp(-X);
rhob = 0*X;
boundary = 0*X;
%boundary(Y>0.5*(ay+by)) = 1;
%Left  & Right
%rhob(1,:) = boundary(1,:); rhob(end,:) = boundary(end,:);
%rhob(1,:) = rhobb(1,:); rhob(end,:) = boundary(end,:);
%Down  &  Up
%rhob(2:end-1,1) = boundary(2:end-1,1); rhob(2:end-1,end) = boundary(2:end-1,end);
rhob=rhobb;
Ihx = speye(Nx+1);
Ihy = speye(Ny+1);
Rxr = Ihx(2:end,:); %D-N
Ryr = Ihy(1:end,:); %N-N
%Rxr = Ihx; %N-N
%Ryr = Ihy; %N-N
%Ryr = [speye(Ny) zeros(Ny,1)]; Ry1(1,1)=1; %Periodic

%% Eigenvalues for Hamiltonian
[Sx1,Lx1] =normeig(Rx1*Ahx*Rx1',Rx1*Bhx*Rx1');
[Sy1,Ly1] = normeig(Ry1*Ahy*Ry1',Ry1*Bhy*Ry1');

[Sx2,Lx2] = normeig(Rx2*Ahx*Rx2',Rx2*Bhx*Rx2');
[Sy2,Ly2] = normeig(Ry2*Ahy*Ry2',Ry2*Bhy*Ry2');

[Sxr,Lxr] =normeig(Rxr*Ahx*Rxr',Rxr*Bhx*Rxr');
[Syr,Lyr] = normeig(Ryr*Ahy*Ryr',Ryr*Bhy*Ryr');

%% Eigenvalues for Shur complement
AAx = Jpx'*Bhx*Dhx*Rx1'*(Rx1*Bhx*Rx1')^(-1)*Rx1*Dhx'*Bhx'*Jpx; AAx = (AAx+AAx')*0.5;
AAy = Jpy'*Bhy*Dhy*Ry2'*(Ry2*Bhy*Ry2')^(-1)*Ry2*Dhy'*Bhy'*Jpy; AAy = (AAy+AAy')*0.5;
BBx = Jpx'*Bhx*Rx1'*(Rx1*Bhx*Rx1')^(-1)*Rx1*Bhx'*Jpx;
BBy = Jpy'*Bhy*Ry2'*(Ry2*Bhy*Ry2')^(-1)*Ry2*Bhy'*Jpy;

[SSx,LLx] =normeig(AAx,BBx);
[SSy,LLy] = normeig(AAy,BBy);

Iphx = speye(size(SSx,1));
Iphy = speye(size(SSy,1));
DD = kron(Iphy,LLx)+kron(LLy,Iphx);
DDinv=full(reshape(diag(DD),size(SSx,1),size(SSy,1))).^(-1);
[i1,i2] = find(DDinv==max(max(abs(DDinv))));
DDinv(i1,i2)=0; %Turn off constant pressure mode

%% Time parameters
c = 2e-2;
dx = min(min(abs(diff(x))),min(abs(diff(y))));
CFL = 0.1;
dt = CFL*dx/c;
tmax = 100;
nsteps = ceil(tmax/dt); 
dt = tmax/nsteps;

%% Time stepping
U10 = U1b; U20 = U2b;
rho0 = rhobb;
%rho0(Y>0.5*(ay+by)) = 1;

Ihx1 = speye(size(Rx1,1));
Ihx2 = speye(size(Rx2,1));
Ihy1 = speye(size(Ry1,1));
Ihy2 = speye(size(Ry2,1));
Ihxr = speye(size(Rxr,1));
Ihyr = speye(size(Ryr,1));

for i=1:nsteps
    %Display
    if mod(i,1000)==0
       disp(sprintf('Iteration %d of %d in realization %d',i,nsteps,n)); 
       %mesh(X,Y,rho1); pause(0.1);
       %mesh(X,Y,U22); pause(0.1);
%        quiver(X,Y,U12,U22); xlabel('x'); ylabel('y'); hold on
%        contour(X,Y,rho1,50); colorbar; hold off ;
%        title(sprintf('N = %d, t = %f', N, i*dt));
%        pause(0.1)
    end
    if i==1 %BDF1/EXT1
        [beta,alpha] = timecoeffs(i);
        Dr = beta(1)*kron(Ihyr,Ihxr)+(nu/Pr)*dt*(kron(Lyr,Ihxr)+kron(Ihyr,Lxr));
        Drinv = full(reshape(diag(Dr),size(Ihxr,1),size(Ihyr,1))).^(-1);
        brtilde = -beta(2)*timestep_op(rho0,Bhx,Bhy); %Add time step from  LHS
        brtilde = brtilde - alpha(1)*dt*advect_op(0,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection
        brtilde = Rxr*((brtilde - inhomogen_op(beta(1),nu,dt,rhob,Ahx,Ahy,Bhx,Bhy))*Ryr'); %Inhomogeneity
        rho1star = Sxr*(Drinv.*(Sxr'*brtilde*Syr))*Syr';
        
        D1 = beta(1)*kron(Ihy1,Ihx1)+nu*dt*(kron(Ly1,Ihx1)+kron(Ihy1,Lx1));
        D1inv = full(reshape(diag(D1),size(Ihx1,1),size(Ihy1,1))).^(-1);
        b1tilde = -beta(2)*timestep_op(U1b,Bhx,Bhy); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(1,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection
        b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        dU1star = Sx1*(D1inv.*(Sx1'*b1tilde*Sy1))*Sy1';
        
        D2 = beta(1)*kron(Ihy2,Ihx2)+nu*dt*(kron(Ly2,Ihx2)+kron(Ihy2,Lx2));
        D2inv = full(reshape(diag(D2),size(Ihx2,1),size(Ihy2,1))).^(-1);
        b2tilde = -beta(2)*timestep_op(U2b,Bhx,Bhy); %Add time step from  LHS;
        b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection;
        b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho0,rho_0,Bhx,Bhy,g); %Add buoyancy
        b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
        dU2star = Sx2*(D2inv.*(Sx2'*b2tilde*Sy2))*Sy2';
       
        [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);

       U11 = Rx1'*U1*Ry1+U1b;
       U21 = Rx2'*U2*Ry2+U2b;
       rho1 = Rxr'*rho1star*Ryr + rhob;
       p1=dp;
       
    elseif i==2 %BDF2/EXT2
        [beta,alpha] = timecoeffs(i);
        [betastar,alphastar] = timecoeffs(i-1);
        
        Dr = beta(1)*kron(Ihyr,Ihxr)+(nu/Pr)*dt*(kron(Lyr,Ihxr)+kron(Ihyr,Lxr));
        Drinv = full(reshape(diag(Dr),size(Ihxr,1),size(Ihyr,1))).^(-1);
        brtilde = -beta(2)*timestep_op(rho1,Bhx,Bhy) - beta(3)*timestep_op(rho0,Bhx,Bhy); %Add time step from  LHS
        brtilde = brtilde - alpha(1)*dt*advect_op(0,U11,U21,rho1,Dhx,Dhy,Bmx,Bmy,Jux,Juy) - alpha(2)*dt*advect_op(0,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection
        brtilde = Rxr*((brtilde - inhomogen_op(beta(1),nu,dt,rhob,Ahx,Ahy,Bhx,Bhy))*Ryr'); %Inhomogeneity
        rho2star = Sxr*(Drinv.*(Sxr'*brtilde*Syr))*Syr';
        
        D1 = beta(1)*kron(Ihy1,Ihx1)+nu*dt*(kron(Ly1,Ihx1)+kron(Ihy1,Lx1));
        D1inv = full(reshape(diag(D1),size(Ihx1,1),size(Ihy1,1))).^(-1);
        b1tilde = -beta(2)*timestep_op(U11,Bhx,Bhy)-beta(3)*timestep_op(U1b,Bhx,Bhy); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(1,U11,U21,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(1,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b1tilde = b1tilde + alphastar(1)*dt*pressure_op(1,p1,Dhx,Bhx,Bhy,Jpx,Jpy);
        b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        dU1star = Sx1*(D1inv.*(Sx1'*b1tilde*Sy1))*Sy1';
       
        D2 = beta(1)*kron(Ihy2,Ihx2)+nu*dt*(kron(Ly2,Ihx2)+kron(Ihy2,Lx2));
        D2inv = full(reshape(diag(D2),size(Ihx2,1),size(Ihy2,1))).^(-1);
        b2tilde = -beta(2)*timestep_op(U21,Bhx,Bhy)-beta(3)*timestep_op(U2b,Bhx,Bhy); %Add time step from  LHS;
        b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U11,U21,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(2,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b2tilde = b2tilde + alphastar(1)*dt*pressure_op(2,p1,Dhy,Bhx,Bhy,Jpx,Jpy);
        b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho1,rho_0,Bhx,Bhy,g) - dt*alpha(2)*buoyancy_op(rho0,rho_0,Bhx,Bhy,g); %Add buoyancy
        b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
        dU2star = Sx2*(D2inv.*(Sx2'*b2tilde*Sy2))*Sy2';
       
       [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);
       

       U12 = Rx1'*U1*Ry1+U1b;
       U22 = Rx2'*U2*Ry2+U2b;
       rho2 = Rxr'*rho2star*Ryr+rhob;
       p2 = alphastar(1)*p1+dp;

        
    else %BDF3/EXT3
        [beta,alpha] = timecoeffs(3);
        [betastar,alphastar] = timecoeffs(2);
        
        Dr = beta(1)*kron(Ihyr,Ihxr)+(nu/Pr)*dt*(kron(Lyr,Ihxr)+kron(Ihyr,Lxr));
        Drinv = full(reshape(diag(Dr),size(Ihxr,1),size(Ihyr,1))).^(-1);
        brtilde = -beta(2)*timestep_op(rho2,Bhx,Bhy) - beta(3)*timestep_op(rho1,Bhx,Bhy) - beta(4)*timestep_op(rho0,Bhx,Bhy); %Add time step from  LHS
        brtilde = brtilde - alpha(1)*dt*advect_op(0,U12,U22,rho2,Dhx,Dhy,Bmx,Bmy,Jux,Juy) - alpha(2)*dt*advect_op(0,U11,U21,rho1,Dhx,Dhy,Bmx,Bmy,Jux,Juy) - alpha(3)*dt*advect_op(0,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection
        brtilde = Rxr*((brtilde - inhomogen_op(beta(1),nu,dt,rhob,Ahx,Ahy,Bhx,Bhy))*Ryr'); %Inhomogeneity
        rho2star = Sxr*(Drinv.*(Sxr'*brtilde*Syr))*Syr';
        
        D1 = beta(1)*kron(Ihy1,Ihx1)+nu*dt*(kron(Ly1,Ihx1)+kron(Ihy1,Lx1));
        D1inv = full(reshape(diag(D1),size(Ihx1,1),size(Ihy1,1))).^(-1);
        b1tilde = -beta(2)*timestep_op(U12,Bhx,Bhy)-beta(3)*timestep_op(U11,Bhx,Bhy)-beta(4)*timestep_op(U10,Bhx,Bhy); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(1,U12,U22,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(1,U11,U21,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)- alpha(3)*dt*advect_op(1,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b1tilde = b1tilde + alphastar(1)*dt*pressure_op(1,p2,Dhx,Bhx,Bhy,Jpx,Jpy)+alphastar(2)*dt*pressure_op(1,p1,Dhx,Bhx,Bhy,Jpx,Jpy);
        b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        dU1star = Sx1*(D1inv.*(Sx1'*b1tilde*Sy1))*Sy1';
        

        D2 = beta(1)*kron(Ihy2,Ihx2)+nu*dt*(kron(Ly2,Ihx2)+kron(Ihy2,Lx2));
        D2inv = full(reshape(diag(D2),size(Ihx2,1),size(Ihy2,1))).^(-1);
        b2tilde = -beta(2)*timestep_op(U22,Bhx,Bhy)-beta(3)*timestep_op(U21,Bhx,Bhy)-beta(4)*timestep_op(U20,Bhx,Bhy); %Add time step from  LHS;
        b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U12,U22,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(2,U11,U21,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)- alpha(3)*dt*advect_op(2,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b2tilde = b2tilde + alphastar(1)*dt*pressure_op(2,p2,Dhy,Bhx,Bhy,Jpx,Jpy)+ alphastar(2)*dt*pressure_op(2,p1,Dhy,Bhx,Bhy,Jpx,Jpy);
        b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho2,rho_0,Bhx,Bhy,g) - dt*alpha(2)*buoyancy_op(rho1,rho_0,Bhx,Bhy,g) - dt*alpha(3)*buoyancy_op(rho0,rho_0,Bhx,Bhy,g); %Add buoyancy
        b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
        dU2star = Sx2*(D2inv.*(Sx2'*b2tilde*Sy2))*Sy2';
        
       [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);
       

       
       U10 = U11; U20 = U21;
       U11  = U12; U21 = U22;
       U12 =Rx1'*U1*Ry1+U1b; U22 = Rx2'*U2*Ry2+U2b;
       
       rho0 = rho1; rho1 = rho2;
       rho2 = Rxr'*rho2star*Ryr+rhob;
       
       pn = alphastar(1)*p2 + alphastar(2)*p1 + dp;
       p1 = p2;
       p2 = pn;
  
    end
end

comp_time = comp_time + toc; comp_time
end

comp_time = comp_time/ntimes; comp_time
% quiver(X,Y,U12,U22); xlabel('x'); ylabel('y'); hold on
% contour(X,Y,rho2,50); colorbar; hold off ;
% title(sprintf('N = %d, t = %f', N, i*dt));
        