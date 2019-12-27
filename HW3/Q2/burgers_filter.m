%
%  Solve Burgers with pseudo-viscosity
%

format compact; format shorte; lw='linewidth';

%if mod(N,2)==0; N=N+1; end;  %%% Enforce N to be odd


c=1; % Speed

h=2*pi/N; x=h*[0:N-1]';



M=N+80;
M=ceil((3*N+1)/2);
C=ones(M,1);         % Constant advecting field

ue=sin(.5*x); ue=ue.^20;
ue=sin(x);

%% SETUP TIMESTEPPER
Tfinal = 6;
CFL=1.00; dxmin=min(abs(diff(x))); dt=CFL*dxmin/c;
Nsteps = ceil(Tfinal/dt);          dt=Tfinal/Nsteps;

iostep = floor(Nsteps/80);

uh=rfft(ue);   %% TRANSFORM DATA INTO WAVE-SPACE
uh=0*uh; uh(3)=1; uh0=uh;

erm=0; iswo=0; time=0; k=0; nsweep = 1;

pad=zeros(M-N,1);

%% Viscous Damping Term --- Diagonal in wave space.
nu=.01; smooth=x; for k=1:N; kw=floor(k/2); smooth(k)=1./(1+nu*dt*kw.^2); end;

for isweep=1:nsweep;
for istep =1:Nsteps; k=k+1;

  u0=uh;       C=irfft([u0; pad]); k1= -dt*advect_op(u0,C,M); 
  u1=u0+.5*k1; C=irfft([u1; pad]); k2= -dt*advect_op(u1,C,M);
  u2=u0+.5*k2; C=irfft([u2; pad]); k3= -dt*advect_op(u2,C,M);
  u3=u0+k3;    C=irfft([u3; pad]); k4= -dt*advect_op(u3,C,M);
  uh=u0 + (k1 + 2*(k2+k3) + k4) / 6; 
  uh=smooth.*uh;                 %%% VISCOUS FILTER
  time = time + dt;

  if mod(istep,iostep)==0;
     u=irfft(uh);
     hold off; plot(x,u,'g-',lw,2)
     mytitle = strcat('1D Spectral Advection: N=',int2str(N));
     title(mytitle,'fontsize',18); drawnow
     pause(.01)
  end;

end;
end;

