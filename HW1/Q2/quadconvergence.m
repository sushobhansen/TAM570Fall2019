maxTquads = zeros(size(nquads,1),1);
for i=1:size(nquads,1)
    [Tr0full,Tpr0full,r,A] = Tsolver(roots,nquads(i),range,L);
    maxTquads(i) = max(abs(sum(Tr0full)));
end
B = -A.*sinh(roots*L)./cosh(roots*L);

%Show convergence of quadrature
figure;
semilogx(nquads,maxTquads);
xlabel('Number of Panels'); 
ylabel('$$\max_{r} |T(r,z=0)|$$','Interpreter','latex');
title('Quadrature Convergence');