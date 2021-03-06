N = 3;
xn = 2*pi*[0:N-1]'/N;
M = 10*N; 
xm = 2*pi*[0:M-1]'/M;

un = [1,4,0]';

um = interp_f(M,un);
plot(xm,um,'r-',xn,un,'bo','linewidth',1.4);
hold on;
plot(xm,0*xm,'k','linewidth',2.0);