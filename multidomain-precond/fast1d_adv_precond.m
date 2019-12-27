format compact; hold off; lw='linewidth';
close all; clear all; clc;
%
% PERIODIC SEM for u_t + c u_x = 0
%
% Set E and N, then run.
%
Nmax=400;E=1;
kappa = zeros(size([1:Nmax]',1),1);
c=1.0;  %% ADVECTING SPEED

ax=-1.0; bx=1.0; lx=bx-ax;

le=lx/E;   xe=ax+(bx-ax)*[0:E]'/E;   % Local (uniform) element length

for N=1:Nmax
N
    [Ah,Bh,Ch,Dh,z,w]=semhat(N);
Q=semq(E,N);  %% No periodicity
nb=size(Q,2); N1=N+1;
R=speye(nb); R=R(2:end,:); %R(1,end)=1; 

QRt=Q*R';   %% Combine Q & R'
QQt=Q*Q';

wq=Q'*ones(N1*E,1); wq=1./wq; %% Counting matrix

% xl=zeros(N1,E);
% for e=1:E; 
%    ae=xe(e); be=xe(e+1); xl(:,e)=ae+0.5*(be-ae)*(1+z);
% end;
% xb=wq.*(Q'*reshape(xl,N1*E,1)); x=xb(1:end-1);

Ae=(2/le)*Ah;
Be=(le/2)*Bh;
Ce=c*Ch;
Bl=zeros(N1,E);
for e=1:E; 
   Bl(:,e)=diag(Be);
end;
Bl=reshape(Bl,N1*E,1); Bl=sparse(Bl); Bl=diag(Bl);
B=QRt'*Bl*QRt;
Bi=1./diag(B);

%% Define matrices
Aspectral = R*Ah*R';
[A1h,B1h,C1h,D1h,z1,w1]=semhat(1);
QFEM = semq(N,1);

dx = diff(z);
Ih = speye(N,N);
AL = kron(spdiags(2./dx,0,N,N),A1h);
AFEM =R*QFEM'*AL*QFEM*R';

kappa(N) = cond(AFEM^(-1)*Aspectral);
plot([1:Nmax]',kappa);
title('D-N Condition'); xlabel('N'); ylabel('\kappa');
axis square; grid on;

% %% INITIAL CONDITION
% ue=(sin(pi*x)).^40;
% 
% %% ESTIMATE DT
% dxmin=min(abs(diff(x))); CFL=0.020; dt=CFL*dxmin/c;

%% SETUP TIMESTEPPER
% Tfinal = 1.;
% Nsteps = ceil(Tfinal/dt); dt=Tfinal/Nsteps;
% 
% iostep = floor(Nsteps/40);
% 
% u=ue; time=0; k=0; nsweep = 1;
% for isweep=1:nsweep;
% for istep =1:Nsteps; k=k+1;
% 
%   u0=u;        k1 = -dt*Bi.*convop_1d(Ce,QRt,u0);
%   u1=u0+.5*k1; k2 = -dt*Bi.*convop_1d(Ce,QRt,u1);
%   u2=u0+.5*k2; k3 = -dt*Bi.*convop_1d(Ce,QRt,u2);
%   u3=u0+k3;    k4 = -dt*Bi.*convop_1d(Ce,QRt,u3);
% 
%   u = u0 + (k1 + 2*(k2+k3) + k4) / 6; 
%   time = time + dt;
% 
% %   if mod(istep,iostep)==0;
% %      ub=R'*u;
% %      hold off; plot(xb,ub,'r.-',lw,2);
% %      mytitle=strcat('1D SEM Convection: (E,N)=(',int2str(E),',',int2str(N),')');
% %      title(mytitle,'fontsize',18); xlabel('- x -','fontsize',18); 
% %      ylabel('- u -','fontsize',18); drawnow
% %      pause(.01)
% %   end;
% 
% end;
% end;
% % ub=R'*u;  plot(xb,ub,'r.-',lw,2); hold on
% % ub=R'*ue; plot(xb,ub,'b--',lw,2);
% 
% err(N) = max(max(abs(u-ue)));
end