clear all; close all; clc;

N=20;

[z,w] = zwgll(N);
R = speye(N+1);R=R(2:end-1,:);
Dh = deriv_mat(z);
Bh = diag(w);
Ah = Dh'*Bh*Dh; Ah=0.5*(Ah+Ah');
A=R*Ah*R';
%f=exp(z);
f = z.*sin(pi*z);
b=R*Bh*f;
u=A\b; ub=[0;u;0];

%exact solution
%a=exp(-1);b=exp(1); ue=(a+(b-a)*.5*(1+z))-exp(z);
ue = (2+2*cos(pi*z)+pi*z.*sin(pi*z))/(pi^3);
err=max(abs(ub-ue));
format shorte; disp([N err]);
plot(z,ub,'ro-',z,ue,'b.-');
legend('Numerical','Exact')