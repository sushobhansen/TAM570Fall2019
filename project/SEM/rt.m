clear; close all; clc;
format longe; format compact;


tol=10e-10;
N=16;
MeshData=readmatrix('Mesh_4.dat','NumHeaderLines',0);
[loc_to_glob,CoorPt,NE]=process_mesh(MeshData,N,tol);
Npt=length(CoorPt(:,1));
Mxy=ceil(1.5*(N+1));
%[loc_to_glob_M,CoorPt_M,NE_M]=process_mesh(MeshData,Mxy,tol);
Np=N-2;
%[loc_to_glob_p,CoorPt_p,NE_p]=process_mesh(MeshData,Np,tol);
figure(1)
plot(CoorPt(:,1),CoorPt(:,2),'o')

%%%%Q matrix test
%u=zeros(Npt,1);
%uL=Qop(u,N,NE,loc_to_glob);
%u2=QTop(uL,N,NE,Npt,loc_to_glob);
%uL=Qop(u2,N,NE,loc_to_glob);
%%%%
lx=zeros(NE,1);
ly=zeros(NE,1);
for ie=1:NE
lx(ie) = CoorPt(loc_to_glob(N+1,ie),1)-CoorPt(loc_to_glob(1,ie),1); 
ly(ie) = CoorPt(loc_to_glob(N*(N+1)+1,ie),2)-CoorPt(loc_to_glob(1,ie),2);
end

Nx=N; Ny=N;
Mx = Mxy; My = Mxy;
Npx =Np; Npy = Np;

nu = 1e-1; g = 9.81; rho_0 = 1000; Pr = 1;
preconflag = 1;
%% U-mesh (P N) and P-mesh (P N-2)
[x,wx] = zwgll(Nx);
[y,wy] = zwgll(Ny);
[xm,wxm] = zwgll(Mx);
[ym,wym] = zwgll(My);
[xp,wpx] = zwgll(Npx);
[yp,wpy] = zwgll(Npy);
X=CoorPt(:,1);
Y=CoorPt(:,2);
Jux = interp_mat(xm,x);
Juy = interp_mat(ym,y);
Jpx = interp_mat(x,xp);
Jpy = interp_mat(y,yp);

Bhx=zeros(Nx+1,Nx+1,NE);
Bhy=zeros(Ny+1,Ny+1,NE);
Bmx=zeros(Mx+1,Mx+1,NE);
Bmy=zeros(My+1,My+1,NE);
Dhx=zeros(Nx+1,Nx+1,NE);
Dhy=zeros(Ny+1,Ny+1,NE);
Ahx=Dhx;
Ahy=Dhy;
for ie=1:NE
Bhx(:,:,ie) = (lx(ie)/2)*spdiags(wx,0,Nx+1,Nx+1);
Bhy(:,:,ie) = (ly(ie)/2)*spdiags(wy,0,Ny+1,Ny+1);
Bmx(:,:,ie) = (lx(ie)/2)*spdiags(wxm,0,Mx+1,Mx+1);
Bmy(:,:,ie) = (ly(ie)/2)*spdiags(wym,0,My+1,My+1);
Dhx(:,:,ie) = (2/lx(ie))*deriv_mat(x);
Dhy(:,:,ie) = (2/ly(ie))*deriv_mat(y);
Ahx(:,:,ie) = Dhx(:,:,ie)'*Bhx(:,:,ie)*Dhx(:,:,ie); Ahx(:,:,ie) = 0.5*(Ahx(:,:,ie)+(Ahx(:,:,ie))');
Ahy(:,:,ie) = Dhy(:,:,ie)'*Bhy(:,:,ie)*Dhy(:,:,ie); Ahy(:,:,ie) = 0.5*(Ahy(:,:,ie)+(Ahy(:,:,ie))');
end

%% Boundaries

%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
boundaryData=zeros((N+1)^2,(N+1)^2,NE);
for ie=1:NE
    boundaryData(:,:,ie)= speye((N+1)^2,(N+1)^2);
end
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
U1b =0*X;
%boundary = 0*X;
%Left  & Right
%U1b(1,:) = boundary(1,:); U1b(end,:) = boundary(end,:);
%U1b(1,:) = 1-y.*y; U1b(end,:) = 1-y.*y;
%Down  &  Up
%U1b(:,1) = boundary(:,1); U1b(:,end) = boundary(:,end);

U2b = 0*Y;

%U2b(1,:) =boundary(1,:); U2b(end,:) = boundary(end,:);
%Down  &  Up
%U2b(:,1) = boundary(:,1); U2b(:,end) = boundary(:,end);

%Mask for boundarys 
ML1=createML(N,NE,boundaryData); %N-N
ML2=createML(N,NE,boundaryData); %N-N


%rhobb = exp(-X);
rhob = exp(-X);

MLr=createML(N,NE,boundaryData); %N-N


%% Time parameters
c = 2e-1;
dx = min(min(abs(diff(x))),min(abs(diff(y))))/2;
CFL = 0.1;
dt = CFL*dx/c;
tmax = 50;
nsteps = ceil(tmax/dt); 
dt = tmax/nsteps;

%% Time stepping
%Initial conditioon
%Gaussian
sigma = 0.1;
U10=exp(-(X.^2+Y.^2)/(2*sigma^2));
U20 = U10;
rho0 = rhob;

for i=1:nsteps
    %Display
    if mod(i,10)==0
       disp(sprintf('Iteration %d of %d',i,nsteps)); 
       %plot3(X,Y,U12,'.')
       %scatter(X,Y,10,U12)
       scatter3(X,Y,U22,10,U22)
       %mesh(X,Y,rho1); pause(0.1);
       %mesh(X,Y,U22); pause(0.1);
       %quiver(X,Y,U12,U22); xlabel('x'); ylabel('y'); hold on
       %contour(X,Y,rho1,50); hold off ;
       pause(0.1)
    end
    if i==1 %BDF1/EXT1
        [beta,alpha] = timecoeffs(i);
        
        rho0L=Qop(rho0,N,NE,loc_to_glob);
        U10L=Qop(U10,N,NE,loc_to_glob);       
        U20L=Qop(U20,N,NE,loc_to_glob);
        
        C10=Maskop(U10L,ML1,NE);
        C20=Maskop(U20L,ML2,NE);

        brtilde = -beta(2)*timestep_semop(rho0L,Bhx,Bhy,NE,N); %Add time step from  LHS
        rhomask0=Maskop(rho0L,MLr,NE);
        brtilde = brtilde - alpha(1)*dt*advect_op(C10,C20,rhomask0,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        %brtilde = Rxr*((brtilde - inhomogen_op(beta(1),nu,dt,rhob,Ahx,Ahy,Bhx,Bhy))*Ryr'); %Inhomogeneity
        brtilde = QTop(Maskop(brtilde,MLr,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        rho1star = solve(preconflag,brtilde,beta(1),Bhx,Bhy,nu/Pr,dt,Ahx,Ahy,NE,N,loc_to_glob,MLr,Npt); 
        

        b1tilde = -beta(2)*timestep_semop(U10L,Bhx,Bhy,NE,N); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(C10,C20,C10,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        b1tilde = QTop(Maskop(b1tilde,ML1,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        dU1star = solve(preconflag,b1tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML1,Npt); 


        b2tilde = -beta(2)*timestep_semop(U20L,Bhx,Bhy,NE,N); %Add time step from  LHS
        b2tilde = b2tilde - alpha(1)*dt*advect_op(C10,C20,C20,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        %b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho0L,rho_0,Bhx,Bhy,g,NE,N); %Add buoyancy
        b2tilde = QTop(Maskop(b2tilde,ML2,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        dU2star = solve(preconflag,b2tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML2,Npt); 
        
       
        %[U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);
        %[U1,U2,dp] = incompress_solve(preconflag,beta(1),dt,AAx,AAy,BBx,BBy,dU1star,dU2star,U1b,U2b,Rx1,Ry1,Rx2,Ry2,Dhx,Dhy,Bhx,Bhy,Jpx,Jpy);
        % [U1,U2,dp] = incompress_semsolve(preconflag,dU1star,dU2star,beta(1),Bhx,Bhy,Dhx,Dhy,nu,dt,Jpx,Jpy,NE,N,Np,loc_to_glob,ML1,ML2,Npt);
       U11 = dU1star;
       U21 = dU2star;
       %rho1 = rho1star + rhob;
       rho1 = rho1star;
       %p1=dp;
       
    elseif i==2 %BDF2/EXT2
        [beta,alpha] = timecoeffs(i);
        [betastar,alphastar] = timecoeffs(i-1);
        
        rho0L=Qop(rho0,N,NE,loc_to_glob);
        rho1L=Qop(rho1,N,NE,loc_to_glob);
        U10L=Qop(U10,N,NE,loc_to_glob);
        U11L=Qop(U11,N,NE,loc_to_glob);
        U20L=Qop(U20,N,NE,loc_to_glob);
        U21L=Qop(U21,N,NE,loc_to_glob);
        
        C11=Maskop(U11L,ML1,NE);
        C21=Maskop(U21L,ML2,NE);
        
        brtilde = -beta(2)*timestep_semop(rho1L,Bhx,Bhy,NE,N)-beta(3)*timestep_semop(rho0L,Bhx,Bhy,NE,N); %Add time step from  LHS
        rhomask0=Maskop(rho0L,MLr,NE);
        rhomask1=Maskop(rho1L,MLr,NE);
        brtilde = brtilde - alpha(1)*dt*advect_op(C11,C21,rhomask1,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(2)*dt*advect_op(C10,C20,rhomask0,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        brtilde = QTop(Maskop(brtilde,MLr,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        rho1star = solve(preconflag,brtilde,beta(1),Bhx,Bhy,nu/Pr,dt,Ahx,Ahy,NE,N,loc_to_glob,MLr,Npt); 
        

        b1tilde = -beta(2)*timestep_semop(U11L,Bhx,Bhy,NE,N)-beta(3)*timestep_semop(U10L,Bhx,Bhy,NE,N); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(C11,C21,C11,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(2)*dt*advect_op(C10,C20,C10,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        b1tilde = QTop(Maskop(b1tilde,ML1,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        dU1star = solve(preconflag,b1tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML1,Npt); 


        b2tilde = -beta(2)*timestep_semop(U21L,Bhx,Bhy,NE,N)-beta(3)*timestep_semop(U20L,Bhx,Bhy,NE,N); %Add time step from  LHS
        b2tilde = b2tilde - alpha(1)*dt*advect_op(C11,C21,C21,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(2)*dt*advect_op(C10,C20,C20,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        %b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho0L,rho_0,Bhx,Bhy,g,NE,N); %Add buoyancy
        b2tilde = QTop(Maskop(b2tilde,ML2,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        dU2star = solve(preconflag,b2tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML2,Npt); 
        
        %%%
        %b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho1,rho_0,Bhx,Bhy,g) - dt*alpha(2)*buoyancy_op(rho0,rho_0,Bhx,Bhy,g); %Add buoyancy
        %%%
       %[U1,U2,dp] = incompress_solve(preconflag,beta(1),dt,AAx,AAy,BBx,BBy,dU1star,dU2star,U1b,U2b,Rx1,Ry1,Rx2,Ry2,Dhx,Dhy,Bhx,Bhy,Jpx,Jpy);

       U12 = dU1star;
       U22 = dU2star;
       rho2 = rho1star;
       %p2 = alphastar(1)*p1+dp;

        
    else %BDF3/EXT3
        [beta,alpha] = timecoeffs(3);
        [betastar,alphastar] = timecoeffs(2);
        
        rho0L=Qop(rho0,N,NE,loc_to_glob);
        rho1L=Qop(rho1,N,NE,loc_to_glob);
        rho2L=Qop(rho2,N,NE,loc_to_glob);
        U10L=Qop(U10,N,NE,loc_to_glob);
        U11L=Qop(U11,N,NE,loc_to_glob);
        U12L=Qop(U12,N,NE,loc_to_glob);
        U20L=Qop(U20,N,NE,loc_to_glob);
        U21L=Qop(U21,N,NE,loc_to_glob);
        U22L=Qop(U22,N,NE,loc_to_glob);
        
        C10=Maskop(U10L,ML1,NE);
        C11=Maskop(U11L,ML2,NE);
        C12=Maskop(U12L,ML2,NE);
        C20=Maskop(U20L,ML1,NE);
        C21=Maskop(U21L,ML2,NE);
        C22=Maskop(U22L,ML2,NE);
        
        brtilde = -beta(2)*timestep_semop(rho2L,Bhx,Bhy,NE,N)-beta(3)*timestep_semop(rho1L,Bhx,Bhy,NE,N)-beta(4)*timestep_semop(rho0L,Bhx,Bhy,NE,N); %Add time step from  LHS
        rhomask0=Maskop(rho0L,MLr,NE);
        rhomask1=Maskop(rho1L,MLr,NE);
        rhomask2=Maskop(rho2L,MLr,NE);
        %brtilde = brtilde - alpha(1)*dt*advect_op(C12,C22,rhomask2,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(2)*dt*advect_op(C11,C21,rhomask1,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(3)*dt*advect_op(C10,C20,rhomask0,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        brtilde = QTop(Maskop(brtilde,MLr,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        rho1star = solve(preconflag,brtilde,beta(1),Bhx,Bhy,nu/Pr,dt,Ahx,Ahy,NE,N,loc_to_glob,MLr,Npt); 
        

        b1tilde = -beta(2)*timestep_semop(U12L,Bhx,Bhy,NE,N)-beta(3)*timestep_semop(U11L,Bhx,Bhy,NE,N)-beta(4)*timestep_semop(U10L,Bhx,Bhy,NE,N); %Add time step from  LHS
        %b1tilde = b1tilde - alpha(1)*dt*advect_op(C12,C22,C12,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(2)*dt*advect_op(C11,C21,C11,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(3)*dt*advect_op(C10,C20,C10,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        b1tilde = QTop(Maskop(b1tilde,ML1,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        dU1star = solve(preconflag,b1tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML1,Npt); 


        b2tilde = -beta(2)*timestep_semop(U22L,Bhx,Bhy,NE,N)-beta(3)*timestep_semop(U21L,Bhx,Bhy,NE,N)-beta(4)*timestep_semop(U20L,Bhx,Bhy,NE,N); %Add time step from  LHS
        b2tilde = b2tilde - alpha(1)*dt*advect_op(C12,C22,C22,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(2)*dt*advect_op(C11,C21,C21,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)- alpha(3)*dt*advect_op(C10,C20,C20,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N); %Add convection
        %b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        %b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho0L,rho_0,Bhx,Bhy,g,NE,N); %Add buoyancy
        b2tilde = QTop(Maskop(b2tilde,ML2,NE),N,NE,Npt,loc_to_glob); %Global restricted btilde
        dU2star = solve(preconflag,b2tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML2,Npt); 
        
      

%         b2tilde = -beta(2)*timestep_op(U22,Bhx,Bhy)-beta(3)*timestep_op(U21,Bhx,Bhy)-beta(4)*timestep_op(U20,Bhx,Bhy); %Add time step from  LHS;
%         b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U12,U22,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(2,U11,U21,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy)- alpha(3)*dt*advect_op(2,U10,U20,rho0,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
%         b2tilde = b2tilde + alphastar(1)*dt*pressure_op(2,p2,Dhy,Bhx,Bhy,Jpx,Jpy)+ alphastar(2)*dt*pressure_op(2,p1,Dhy,Bhx,Bhy,Jpx,Jpy);
%         b2tilde = b2tilde - dt*alpha(1)*buoyancy_op(rho2,rho_0,Bhx,Bhy,g) - dt*alpha(2)*buoyancy_op(rho1,rho_0,Bhx,Bhy,g) - dt*alpha(3)*buoyancy_op(rho0,rho_0,Bhx,Bhy,g); %Add buoyancy
%         b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
%         dU2star = solve(preconflag,b2tilde,beta(1),Bhx,Bhy,nu,dt,Ahx,Ahy,Rx2,Ry2);
        
       %[U1fdm,U2fdm,dpfdm] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);
       %[U1,U2,dp] = incompress_solve(preconflag,beta(1),dt,AAx,AAy,BBx,BBy,dU1star,dU2star,U1b,U2b,Rx1,Ry1,Rx2,Ry2,Dhx,Dhy,Bhx,Bhy,Jpx,Jpy);

       
       U10 = U11; U20 = U21;
       U11  = U12; U21 = U22;
       U12 =dU1star; U22 = dU2star;
       
       rho0 = rho1; rho1 = rho2;
       rho2 = rho1star;
       
       %pn = alphastar(1)*p2 + alphastar(2)*p1 + dp;
       %p1 = p2;
       %p2 = pn;
  
    end
end

% streamslice(X',Y',U12',U22');
% xlim([ax,bx]); ylim([ay,by]);
% xlabel('x'); ylabel('y');
% title(sprintf('Streamlines for N=%d, Re=%d, Time=%f',N,Re,tmax));
% axis square;
% quiver(X,Y,U12,U22); hold on; xlabel('x'); ylabel('y');
% contourf(X,Y,vort(U12,U22,Dhx,Dhy),10); title('Vorticity and velocity contours');
% figure;
% contourf(X,Y,Jpx*p2*Jpy');colorbar
% xlabel('x'); ylabel('y'); title('Pressure');
        