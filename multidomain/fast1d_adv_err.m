format compact; lw='linewidth';
 clear all; clc;
%
% PERIODIC SEM for u_t + c u_x = 0
%
% Set E and N, then run.
%
Nmax=40;E=5;
S = zeros(size([1:Nmax]',1),1);
err = zeros(size([1:Nmax]',1),1);
c=1.0;  %% ADVECTING SPEED

ax=0; bx=1.0; lx=bx-ax;

le=lx/E;   xe=ax+(bx-ax)*[0:E]'/E;   % Local (uniform) element length

for N=1:Nmax
N
    [Ah,Bh,Ch,Dh,z,w]=semhat(N);
Q=semq(E,N);  %% No periodicity
nb=size(Q,2); N1=N+1;
R=speye(nb); R=R(1:end-1,:); R(1,end)=1; 

QRt=Q*R';   %% Combine Q & R'
QQt=Q*Q';

wq=Q'*ones(N1*E,1); wq=1./wq; %% Counting matrix

xl=zeros(N1,E);
for e=1:E; 
   ae=xe(e); be=xe(e+1); xl(:,e)=ae+0.5*(be-ae)*(1+z);
end;
xb=wq.*(Q'*reshape(xl,N1*E,1)); x=xb(1:end-1);

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


%% INITIAL CONDITION
ue=(sin(pi*x)).^40;

%% ESTIMATE DT
dxmin=min(abs(diff(x))); CFL=0.020; dt=CFL*dxmin/c;

%% GET OPERATOR BiC
% Ih = eye(E*N);
% L = zeros(E*N,E*N);
% for j=1:E*N
%    L(:,j) = -Bi.*convop_1d(Ce,QRt,Ih(:,j)); 
% end
% 
% D = eig(L);
% Dmax = max(abs(D));
% S(N) = (c/dxmin)/Dmax;

% SETUP TIMESTEPPER
Tfinal = 1.;
Nsteps = ceil(Tfinal/dt); dt=Tfinal/Nsteps;

iostep = floor(Nsteps/40);

u=ue; time=0; k=0; nsweep = 1;
for isweep=1:nsweep;
for istep =1:Nsteps; k=k+1;

  u0=u;        k1 = -dt*Bi.*convop_1d(Ce,QRt,u0);
  u1=u0+.5*k1; k2 = -dt*Bi.*convop_1d(Ce,QRt,u1);
  u2=u0+.5*k2; k3 = -dt*Bi.*convop_1d(Ce,QRt,u2);
  u3=u0+k3;    k4 = -dt*Bi.*convop_1d(Ce,QRt,u3);

  u = u0 + (k1 + 2*(k2+k3) + k4) / 6; 
  time = time + dt;

%   if mod(istep,iostep)==0;
%      ub=R'*u;
%      hold off; plot(xb,ub,'r.-',lw,2);
%      mytitle=strcat('1D SEM Convection: (E,N)=(',int2str(E),',',int2str(N),')');
%      title(mytitle,'fontsize',18); xlabel('- x -','fontsize',18); 
%      ylabel('- u -','fontsize',18); drawnow
%      pause(.01)
%   end;

end;
end;
% ub=R'*u;  plot(xb,ub,'r.-',lw,2); hold on
% ub=R'*ue; plot(xb,ub,'b--',lw,2);

err(N) = max(max(abs(u-ue)));
end;

semilogy(1:Nmax,err,'.','LineWidth',3); 
xlabel('N'); ylabel('Error');