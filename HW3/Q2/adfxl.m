clear all; close all; clc;

Nx = 31 ;
lx = 2*pi;
hx = lx/Nx; x = hx*[0:Nx-1]';
Mx = ceil(1/.5*(Nx+1));
Df = dhatf(Nx);
Dff = Df^2;
nu = 0.1;
Bh = pi*speye(Nx); Bh(1,1) = 2*pi;

u0 = zeros(size(x));
u0(x>=1 & x<=2) = 1;
uxm = interp_Df(Mx,Df,u0);
cx = 1+0*x;
cm = interp_f(Mx,cx);
wm = cm.*uxm;

P = [speye(Nx); zeros(Mx-Nx,Nx)];

cfl = 0.1; c=1;
dt = cfl*hx/c;
tmax = 2*pi;
nsteps = ceil(tmax/dt); dt = tmax/nsteps;

u = u0;
for i=1:nsteps
    u = rfft(u);
    u = (speye(Nx) - Bh^(-1)*dt*nu*Dff)\(u - dt*P'*rfft(wm));
    u = irfft(u);
    wm = cm.*interp_Df(Mx,Df,u);
    if(mod(i,50)==0)
    plot(x,u);
    title(sprintf('Time = %f',i*dt));
    pause();
    end
end