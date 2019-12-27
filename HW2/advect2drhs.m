function V = advect2drhs(Bh,Bm,Cxm,Cym,Dtilde,Jh,R,U)
%Returns -B^(-1)*C*u with appropriate boundary conditions applied as
%defined by R
m = size(Cxm,1);
n = size(U,1);
Bh = R*Bh*R';
U = R*U*R';
V = Cxm.*((Dtilde*R')*U*(R*Jh')) + Cym.*((Jh*R')*U*(R*Dtilde'));
V = reshape(full(diag(Bm)),m,m).*V;
V = (R*Jh')*V*(Jh*R'); %CU
V = -(Bh^(-1))*V*(Bh^(-1))';
V = R'*V*R;

end

