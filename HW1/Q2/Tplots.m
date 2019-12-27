Tr0 = zeros(n1,size(r,1));
Tpr0 = zeros(n1,size(r,1));
for j=1:n1
   Tr0(j,:) = sum(Tr0full(1:j,:),1);
   Tpr0(j,:) = sum(Tpr0full(1:j,:),1); 
end

figure; 
for j=1:n1
    plot(r,Tr0(j,:),'DisplayName',sprintf('N = %d',j)); 
    hold on;
end
xlabel('r'); ylabel('$T(r,z=0)$','Interpreter','latex'); title(sprintf('Number of Panels = %d',nquads(end)));
legend;

figure;
for j=1:n1
    plot(r,Tpr0(j,:),'DisplayName',sprintf('N = %d',j));
    hold on;
end
xlabel('r'); ylabel('$T_z(r,z=0)$','Interpreter','latex'); title(sprintf('Number of Panels = %d',nquads(end)));
legend('Location','southwest');

z = linspace(0,L,100);
[R,Z] = meshgrid(r,z);
T = zeros(size(R));

for k=1:n
   T = T + (A(k)*sinh(roots(k)*Z)-A(k)*cosh(roots(k)*Z)).*besselj(0,roots(k)*R);
end
figure;
contourf(R,Z,T,50,'edgecolor','none'); xlabel('r'); ylabel('z'); zlabel('T'); ylim([0,3]); colorbar hot;