clear all;
a = 0.001;

L = 1.0; %size of domain
n = 500;

c = 1.0;

h = L/n;
n1 = n-1;
e = ones(n1,1);


A = -a*spdiags([e -2*e e], -1:1, n1, n1)/(h*h); %A - SPD
C = c*spdiags([-e 0*e e], -1:1, n1, n1)/(2*h); %C - skew-symm
H = A+C; %AD operator

x = [1:n1]*h; 
x = x';

f = 1.0+0*x;
u = H\f;
ub = [0;u;0];
xb = [0;x;L]; %Extend u and x to boundaries

%Analytical solution
eL = exp(-c*L/a);
ex = exp(c*(xb-L)/a);
utilde = (1/c)*(xb-(L*(ex-eL))/(1-eL));

%Numerical, closed form
Peg = c*h/a;
gamma = (2+Peg)/(2-Peg);
unr = (L/c)*(x/L - (1-gamma.^[1:n1]')/(1-gamma^n));
unrb = [0;unr;0];

plot(xb,ub,'b-','DisplayName','Numerical (Eq 31)'); hold on;
plot(xb,utilde,'r*','DisplayName','Analytical (Eq 30)');
plot(xb,unrb,'kd','DisplayName','Numerical (closed form, Eq 32)');
xlabel('x'); ylabel('u'); title(sprintf('Pe = %d, n = %d',c*L/a, n));
legend('Location','south');