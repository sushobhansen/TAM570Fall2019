clear all;
a = [1.0,0.1,0.01,0.001]';

L = 1.0; %size of domain
n = 2.^[1:20]';

c = 1.0;

for i=1:size(a,1)
error = AD1Dsolver(a(i),c,L,n);

loglog(n,error,'DisplayName',sprintf('Pe = %d',c*L/a(i)));
hold on;
end

xlabel('n'); ylabel('scaled error'); title('1D Advection-Diffusion Convergence');
legend;