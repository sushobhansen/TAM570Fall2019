clear; close all; clc;
format longe; format compact;

X1=-0.5; X1=-1; 
X2=6.0; X2=10;
Y1=-0.5; Y1=-1;
Y2=1.5; Y2=1;
lx = X2-X1; ly = Y2-Y1;
Nx = 30; Ny = 15;
Mx = ceil(1.5*(Nx+1)); My = ceil(1.5*(Ny+1));
Npx = Nx-2; Npy = Ny-2;
Re = 200; nu = 1/Re; %nu=0.001;
tol=10e-8;
randflag=0;%random initial condition 
global kovasznaycase;%kovasznay case: plus 1 in velocity during advection collocation.
kovasznaycase=1;

%% U-mesh (P N)
[x,wx] = zwgll(Nx);x=X1+(x+1)*lx/2;
[y,wy] = zwgll(Ny);y=Y1+(y+1)*ly/2;
[xm,wxm] = zwgll(Mx);xm=X1+(xm+1)*lx/2;
[ym,wym] = zwgll(My);ym=Y1+(ym+1)*ly/2;
Bhx = (lx/2)*spdiags(wx,0,Nx+1,Nx+1);
Bhy = (ly/2)*spdiags(wy,0,Ny+1,Ny+1);
Bmx = (lx/2)*spdiags(wxm,0,Mx+1,Mx+1);
Bmy = (ly/2)*spdiags(wym,0,My+1,My+1);
[Xplot,Yplot]=meshgrid(x,y);
[X,Y] = ndgrid(x,y);
Dhx = (2/lx)*deriv_mat(x);
Dhy = (2/ly)*deriv_mat(y);
Ahx = Dhx'*Bhx*Dhx; Ahx = 0.5*(Ahx+Ahx');
Ahy = Dhy'*Bhy*Dhy; Ahy = 0.5*(Ahy+Ahy');
Jux = interp_mat(xm,x);
Juy = interp_mat(ym,y);

%%
Xex=X;
Yex=Y;
lambda2=-4*pi^2/Re;
u1ex=-exp(lambda2*Xex).*cos(2*pi*Yex);
u2ex=(lambda2/(2*pi))*exp(lambda2*Xex).*sin(2*pi*Yex);

%% P-mesh (P N-2)
[xp,wpx] = zwgll(Npx);xp=X1+(xp+1)*lx/2;
[yp,wpy] = zwgll(Npy);yp=Y1+(yp+1)*ly/2;
Jpx = interp_mat(x,xp);
Jpy = interp_mat(y,yp);

%% Boundaries
c = 1.;
U1b = 0*X; 
%boundary1 = -exp(lambda2*(X)).*cos(2*pi*(Y));
boundary1 = 1;

%Left  & Right
U1b(1,:) =boundary1(1,:);  U1b(end,:) =0;
%U1b(1,:) =0; U1b(end,:) =0;

%Down  &  Up
U1b(:,1) = boundary1(:,1); U1b(:,end) = boundary1(:,end);


U2b = 0*X; 
%boundary2=0.5*(lambda2/pi)*exp(lambda2*(X)).*sin(2*pi*(Y));
boundary2=0*Y;
%Left & Right
%U2b(1,:)=(lambda2/(2*pi))*exp(lambda2*(-0.5)).*sin(2*pi*(y+0.5)); U2b(end,:)=0;
U2b(1,:)=boundary2(1,:); U2b(end,:)=0;
%Down & Up
U2b(:,1)=boundary2(:,1); U2b(:,end)=boundary2(:,end);

%U2b=U1b; %U1b=0*X;


Ihx = speye(Nx+1);
Ihy = speye(Ny+1);
%Rx1 = Ihx(2:end-1,:);%D-D
Rx1 = Ihx(2:end,:);%Left inhomo D, right homo N
Ry1 = Ihy(2:end-1,:);%D-D
%Ry1 = Ihy(1:end,:);%Upper N Lower N
%Rx1 = [speye(Nx); zeros(1,Nx)]; Rx1(end,1)=1;Rx1=Rx1'; %Periodic
%Ry1 = [speye(Ny); zeros(1,Ny)]; Ry1(end,1)=1;Ry1=Ry1'; %Periodic

%Rx2 = Ihx(2:end-1,:); %D-D
Rx2 = Ihx(2:end,:); %Left inhomo D, right homo N
Ry2 = Ihy(2:end-1,:);%D-D
%Ry2 = Ihy(1:end,:);%Upper N Lower N
%Rx2 = [speye(Nx); zeros(1,Nx)]; Rx2(end,1)=1;Rx2=Rx2'; %Periodic
%Ry2 = [speye(Ny); zeros(1,Ny)]; Ry2(end,1)=1; Ry2=Ry2';%Periodic

%% Eigenvalues for Hamiltonian
Bx1=Rx1*Bhx*Rx1';
By1=Ry1*Bhy*Ry1';
Bx2=Rx2*Bhx*Rx2';
By2=Ry2*Bhy*Ry2';
Ahrx1 = (Rx1*Dhx'*Rx1')*(Rx1*Bhx*Rx1')*(Rx1*Dhx*Rx1'); Ahrx1  = 0.5*(Ahrx1+Ahrx1');
Ahry1 = (Ry1*Dhy'*Ry1')*(Ry1*Bhy*Ry1')*(Ry1*Dhy*Ry1'); Ahry1  = 0.5*(Ahry1+Ahry1');
Ahrx2 = (Rx2*Dhx'*Rx2')*(Rx2*Bhx*Rx2')*(Rx2*Dhx*Rx2'); Ahrx2  = 0.5*(Ahrx2+Ahrx2');
Ahry2 = (Ry2*Dhy'*Ry2')*(Ry2*Bhy*Ry2')*(Ry2*Dhy*Ry2'); Ahry2  = 0.5*(Ahry2+Ahry2');
%[Sx1,Lx1] = eig(full(Rx1*Ahx*Rx1'),full(Rx1*Bhx*Rx1'));
[Sx1,Lx1] =normeig(Rx1*Ahx*Rx1',Bx1);
%[Sy1,Ly1] = eig(full(Ry1*Ahy*Ry1'),full(Ry1*Bhy*Ry1'));
[Sy1,Ly1] = normeig(Ry1*Ahy*Ry1',By1);
Ly1(1,1)=abs(Ly1(1,1));
%[Sx2,Lx2] = eig(full(Rx2*Ahx*Rx2'),full(Rx2*Bhx*Rx2'));
[Sx2,Lx2] = normeig(Rx2*Ahx*Rx2',Bx2);

%[Sy2,Ly2] = eig(full(Ry2*Ahy*Ry2'),full(Ry2*Bhy*Ry2'));
[Sy2,Ly2] = normeig(Ry2*Ahy*Ry2',By2);
%Sy2 = bsxfun(@rdivide,Sy2,vecnorm(Sy2'*(Ry2*Bhy*Ry2')*Sy2));

Ly2(1,1)=abs(Ly2(1,1));
%% Eigenvalues for Shur complement
%AAx = Jpx'*Bhx*Dhx*Rx1'*(Rx1*Bhx*Rx1')^(-1)*Rx1*Dhx'*Bhx'*Jpx; AAx = (AAx+AAx')*0.5;
%AAy = Jpy'*Bhy*Dhy*Ry2'*(Ry2*Bhy*Ry2')^(-1)*Ry2*Dhy'*Bhy'*Jpy; AAy = (AAy+AAy')*0.5;
%BBx = Jpx'*Bhx*Rx1'*(Rx1*Bhx*Rx1')^(-1)*Rx1*Bhx'*Jpx;
%BBy = Jpy'*Bhy*Ry2'*(Ry2*Bhy*Ry2')^(-1)*Ry2*Bhy'*Jpy;
AAx = Jpx'*Bhx*Dhx*Rx1'*(Rx1*Bhx*Rx1')^(-1)*Rx1*Dhx'*Bhx'*Jpx; AAx = (AAx+AAx')*0.5;
AAy = Jpy'*Bhy*Dhy*Ry2'*(Ry2*Bhy*Ry2')^(-1)*Ry2*Dhy'*Bhy'*Jpy; AAy = (AAy+AAy')*0.5;
BBx = Jpx'*Bhx*Rx1'*(Rx1*Bhx*Rx1')^(-1)*Rx1*Bhx'*Jpx;
BBy = Jpy'*Bhy*Ry2'*(Ry2*Bhy*Ry2')^(-1)*Ry2*Bhy'*Jpy;

[SSx,LLx] =normeig(AAx,BBx);

%[SSx,LLx] = eig(AAx,BBx);
%normalizer = vecnorm(SSx'*BBx*SSx)';
%for i=1:size(SSx,2)
%    SSx(:,i) = SSx(:,i)./sqrt(SSx(:,i)'*BBx*SSx(:,i));
%end
%SSx = bsxfun(@rdivide,SSx,vecnorm(SSx'*BBx*SSx));
[SSy,LLy] = normeig(AAy,BBy);
%SSy = bsxfun(@rdivide,SSy,vecnorm(SSy'*BBy*SSy));
%for i=1:size(SSy,2)
%    SSy(:,i) = SSy(:,i)./sqrt(SSy(:,i)'*BBy*SSy(:,i));
%end
Iphx = speye(size(SSx,1));
Iphy = speye(size(SSy,1));
DD = kron(Iphy,LLx)+kron(LLy,Iphx);
if kovasznaycase==0
DD = full(reshape(diag(DD),size(SSx,1),size(SSy,1)));
[i1,i2] = find(DD==min(min(DD)));
DDinv = DD.^(-1); DDinv(i1,i2)=0; %Turn off constant pressure mode
elseif kovasznaycase==1
DDinv = DD^(-1);
DDinv=full(reshape(diag(DDinv),size(SSx,1),size(SSy,1)));
else
    disp('ERROR!!! Wrong kovasznaycase');
    return
end
%% Time parameters

dx = min(min(abs(diff(x))),min(abs(diff(y))));
CFL = 0.5;
dt = CFL*dx/c;
tmax = 20;
nsteps = ceil(tmax/dt); 
dt = tmax/nsteps;
err1=zeros(nsteps,1);
err2=err1;
errp=err1;
%% Time stepping
U10 = U1b + cos(0.5*pi*X).*sin(0.5*pi*Y); U20 = U2b - sin(0.5*pi*X).*cos(0.5*pi*Y);
U11 =randn(size(U1b))*randflag;U21=randn(size(U2b))*randflag;
U12 =randn(size(U1b))*randflag+(randflag-1);U22=randn(size(U2b))*randflag+(randflag-1);
p2=0*X(2:end-1,2:end-1)+1;p1=0*X(2:end-1,2:end-1)+0.5;
Ihx1 = speye(size(Rx1,1));
Ihx2 = speye(size(Rx2,1));
Ihy1 = speye(size(Ry1,1));
Ihy2 = speye(size(Ry2,1));
%nsteps=50;
dU2star=0*X(2:end,1:end-1);
U2=dU2star;
for i=1:nsteps
    if mod(i,10)==0
        i
        %contourf(X,Y,Jpx*p2*Jpy');colorbar;hold on
        contourf(X,Y,vort(U12,U22,Dhx,Dhy),30);colorbar;hold on;
        %contourf(X,Y,(U12.^2+U22.^2).^0.5); hold on
        quiver(X,Y,U12,U22); xlabel('x'); ylabel('y');hold off;
        axis equal; pause(0.1);
        %surf(X,Y,(U12.^2+U22.^2).^0.5)
        %surf(X,Y,U22)
        %plot(y(1:end-1),dU2star(1,:)/max(dU2star(1,:)));hold on;
        %plot(y(1:end-1),U2(1,:)/max(U2(1,:)));
        %plot(y(1:end-1),dU2star(1,:));hold on;
        %plot(y(1:end-1),U2(1,:));
        %legend('0star','0','1star','1','2star','2','3star','3')
        %surf(X,Y,Jpx*p2*Jpy');colorbar;
       pause(0.1);
    end
    if i==1
        [beta,alpha] = timecoeffs(i);
        D1 = beta(1)*kron(Ihy1,Ihx1)+nu*dt*(kron(Ly1,Ihx1)+kron(Ihy1,Lx1));
        D1inv=D1^-1;
        D1inv = full(reshape(diag(D1inv),size(Ihx1,1),size(Ihy1,1)));
        b1tilde = -beta(2)*timestep_op(U1b,Bhx,Bhy); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(1,U10,U20,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection
        b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        dU1star = Sx1*(D1inv.*(Sx1'*b1tilde*Sy1))*Sy1';
        
        D2 = beta(1)*kron(Ihy2,Ihx2)+nu*dt*(kron(Ly2,Ihx2)+kron(Ihy2,Lx2));
        D2inv=D2^-1;
        D2inv = full(reshape(diag(D2inv),size(Ihx2,1),size(Ihy2,1)));
        b2tilde = -beta(2)*timestep_op(U2b,Bhx,Bhy); %Add time step from  LHS;
        b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U10,U20,Dhx,Dhy,Bmx,Bmy,Jux,Juy); %Add convection;
        b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
        dU2star = Sx2*(D2inv.*(Sx2'*b2tilde*Sy2))*Sy2';
       if kovasznaycase==2
           disp('kovasznaycase First step p overwirte with analytical soln')
           dU1star=u1ex(2:end,1:end-1);
           dU2star=u2ex(2:end,1:end-1);
       end
       [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);

       U11 = Rx1'*U1*Ry1+U1b;
       U21 = Rx2'*U2*Ry2+U2b;
       if kovasznaycase==2
           disp('kovasznaycase First step overwirte with analytical soln')
           U11=u1ex;
           U21=u2ex;
       end
       p1 = dp;
    elseif i==2
        [beta,alpha] = timecoeffs(i);
        [betastar,alphastar] = timecoeffs(i-1);
        D1 = beta(1)*kron(Ihy1,Ihx1)+nu*dt*(kron(Ly1,Ihx1)+kron(Ihy1,Lx1));
        D1inv=D1^-1;
        D1inv = full(reshape(diag(D1inv),size(Ihx1,1),size(Ihy1,1)));
        b1tilde = -beta(2)*timestep_op(U11,Bhx,Bhy)-beta(3)*timestep_op(U1b,Bhx,Bhy); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(1,U11,U21,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(1,U10,U20,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b1tilde = b1tilde + alphastar(1)*dt*pressure_op(1,p1,Dhx,Bhx,Bhy,Jpx,Jpy);
        b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        dU1star = Sx1*(D1inv.*(Sx1'*b1tilde*Sy1))*Sy1';
       
        D2 = beta(1)*kron(Ihy2,Ihx2)+nu*dt*(kron(Ly2,Ihx2)+kron(Ihy2,Lx2));
        D2inv=D2^-1;
        D2inv = full(reshape(diag(D2inv),size(Ihx2,1),size(Ihy2,1)));
        b2tilde = -beta(2)*timestep_op(U21,Bhx,Bhy)-beta(3)*timestep_op(U2b,Bhx,Bhy); %Add time step from  LHS;
        b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U11,U21,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(2,U10,U20,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b2tilde = b2tilde + alphastar(1)*dt*pressure_op(2,p1,Dhy,Bhx,Bhy,Jpx,Jpy);
        b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
        dU2star = Sx2*(D2inv.*(Sx2'*b2tilde*Sy2))*Sy2';
       if kovasznaycase==2
           disp('kovasznaycase Second step p overwirte with analytical soln')
           dU1star=u1ex(2:end,1:end-1);
           dU2star=u2ex(2:end,1:end-1);
       end
       [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);
       

       U12 = Rx1'*U1*Ry1+U1b;
       U22 = Rx2'*U2*Ry2+U2b;
       if kovasznaycase==2
           disp('kovasznaycase Second step overwirte with analytical soln')
           U12=u1ex;
           U22=u2ex;
       end
       p2 = alphastar(1)*p1+dp;

        
    else
        [beta,alpha] = timecoeffs(3);
        [betastar,alphastar] = timecoeffs(2);
        
        D1 = beta(1)*kron(Ihy1,Ihx1)+nu*dt*(kron(Ly1,Ihx1)+kron(Ihy1,Lx1));
        D1inv=D1^-1;
        D1inv = full(reshape(diag(D1inv),size(Ihx1,1),size(Ihy1,1)));
        b1tilde = -beta(2)*timestep_op(U12,Bhx,Bhy)-beta(3)*timestep_op(U11,Bhx,Bhy)-beta(4)*timestep_op(U10,Bhx,Bhy); %Add time step from  LHS
        b1tilde = b1tilde - alpha(1)*dt*advect_op(1,U12,U22,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(1,U11,U21,Dhx,Dhy,Bmx,Bmy,Jux,Juy)- alpha(3)*dt*advect_op(1,U10,U20,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b1tilde = b1tilde + alphastar(1)*dt*pressure_op(1,p2,Dhx,Bhx,Bhy,Jpx,Jpy)+alphastar(2)*dt*pressure_op(1,p1,Dhx,Bhx,Bhy,Jpx,Jpy);
        b1tilde = Rx1*((b1tilde - inhomogen_op(beta(1),nu,dt,U1b,Ahx,Ahy,Bhx,Bhy))*Ry1'); %Inhomogeneity
        dU1star = Sx1*(D1inv.*(Sx1'*b1tilde*Sy1))*Sy1';
        

        D2 = beta(1)*kron(Ihy2,Ihx2)+nu*dt*(kron(Ly2,Ihx2)+kron(Ihy2,Lx2));
        D2inv=D2^-1;
        D2inv = full(reshape(diag(D2inv),size(Ihx2,1),size(Ihy2,1)));
        b2tilde = -beta(2)*timestep_op(U22,Bhx,Bhy)-beta(3)*timestep_op(U21,Bhx,Bhy)-beta(4)*timestep_op(U20,Bhx,Bhy); %Add time step from  LHS;
        b2tilde = b2tilde - alpha(1)*dt*advect_op(2,U12,U22,Dhx,Dhy,Bmx,Bmy,Jux,Juy)-alpha(2)*dt*advect_op(2,U11,U21,Dhx,Dhy,Bmx,Bmy,Jux,Juy)- alpha(3)*dt*advect_op(2,U10,U20,Dhx,Dhy,Bmx,Bmy,Jux,Juy);
        b2tilde = b2tilde + alphastar(1)*dt*pressure_op(2,p2,Dhy,Bhx,Bhy,Jpx,Jpy)+ alphastar(2)*dt*pressure_op(2,p1,Dhy,Bhx,Bhy,Jpx,Jpy);
        b2tilde = Rx2*((b2tilde - inhomogen_op(beta(1),nu,dt,U2b,Ahx,Ahy,Bhx,Bhy))*Ry2');
        dU2star = Sx2*(D2inv.*(Sx2'*b2tilde*Sy2))*Sy2';
        
       [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta(1),U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2);
       

       
       U10 = U11; U20 = U21;
       U11  = U12; U21 = U22;
       U12 =Rx1'*U1*Ry1+U1b; U22 = Rx2'*U2*Ry2+U2b;
       
       pn = alphastar(1)*p2 + alphastar(2)*p1 + dp;
       p1 = p2;
       p2 = pn;
       
  
    end
    err1(i)=max(max(abs(U12-u1ex)));
    err2(i)=max(max(abs(U22-u2ex)));
    errp(i)=max(max(abs(p2-p1)));
    if err1(i)+err2(i)+errp(i)<tol
        break
    elseif err1(i)+err2(i)>100.
        disp('ERROR!!! Result blows up!!!!');
        mesh(X,Y,Jpx*p2*Jpy')
        return
    end
end
        figure(1)
       contourf(X,Y,vort(U11,U21,Dhx,Dhy),30);
       hold on
       quiver(X,Y,U11,U21); hold on; xlabel('x'); ylabel('y');
        hold off
       figure(2);
       contourf(X,Y,Jpx*p2*Jpy');colorbar
        xlabel('x'); ylabel('y');
        
       figure(3)
       loglog(err1)
       hold on
       loglog(err2)
       loglog(errp)
       hold off
       legend('Relative error of U_1','Relative error of U_2','Relative error of pressure')
       
       figure(4)

       contourf(X,Y,sqrt(u1ex.^2+u2ex.^2)); hold on; xlabel('x'); ylabel('y');colorbar; title('Exact');
       
       figure(5)
       contourf(X,Y,sqrt(U12.^2+U22.^2)); hold on; xlabel('x'); ylabel('y');colorbar; title('Numerical');
       
       figure(6)
       contourf(X,Y,abs(sqrt((U12-u1ex).^2+(U22-u2ex).^2))); hold on; xlabel('x'); ylabel('y');colorbar; title('Error');
       