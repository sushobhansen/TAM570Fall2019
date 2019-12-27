Tmax = -sum(A)*besselj(0,0);
TN = zeros(n,1);
Tsum = 0;

for k=1:n
   Tsum = Tsum + (A(k)*sinh(roots(k)*Z)-A(k)*cosh(roots(k)*Z)).*besselj(0,roots(k)*R);
   TN(k) = max(Tsum,[],'all');
end

figure;
loglog(1:n,abs(TN-Tmax));
xlabel('N'); ylabel('$|T_N - T_{max}|$','Interpreter','latex');