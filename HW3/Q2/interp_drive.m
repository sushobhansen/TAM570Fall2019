

N=3;        xn=2*pi*[0:N-1]'/N;
M=10*(N+5); xm=2*pi*[0:M-1]'/M;

un = [1; 4; 0];

um = interp_f(M,un);

plot(xm,0*xm,'k-',xm,um,'r-',xn,un,'bo','linewidth',1.4)

pause;



N=4;        xn=2*pi*[0:N-1]'/N;
M=10*(N+5); xm=2*pi*[0:M-1]'/M;

un = sin(xn);

um = interp_f(M,un);

plot(xm,0*xm,'k-',xm,um,'r-',xn,un,'bo','linewidth',1.4)

pause;



N=81;       xn=2*pi*[0:N-1]'/N;
M=10*(N+5); xm=2*pi*[0:M-1]'/M;

un = sign(xn-pi);

um = interp_f(M,un);

plot(xm,0*xm,'k-',xm,um,'r-',xn,un,'bo','linewidth',1.4)

pause;

