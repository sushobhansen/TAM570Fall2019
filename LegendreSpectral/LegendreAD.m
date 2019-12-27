clear all; close all; clc;

a = 0.01; c = 1.0; L=1.0;
Pe = c*L/a;
Nmax = 10;
error_sem = zeros(Nmax,1);
error_fd = zeros(Nmax,1);

for i=1:Nmax
    N = 2^i;
    [z,w] = zwgll(N);
    z = (z+1)*0.5;
    w = w*0.5;    
    numr = exp(Pe*(z-1))-exp(-Pe);
    denr = 1-exp(-Pe);
    utilde = (L/c)*(z-(numr/denr));
    %Legendre Spectral Method
    %Restriction matrix
    R = speye(N+1);
    R = R(2:end-1,:);

    Dh = deriv_mat(z);
    Bh = spdiags(w,0,N+1,N+1);
    Ah = a*Dh'*Bh*Dh; Ah = 0.5*(Ah+Ah');
    A = R*Ah*R';
    Ch = c*Bh*Dh;
    C = R*Ch*R';

    f = ones(size(z));
    b = R*Bh*f;

    u = (A+C)\b;
    ub = [0;u;0];
    plot(z,ub,'DisplayName',sprintf('N = %d',2^i)); hold on;
    error_sem(i) = max(abs(ub-utilde));
    
    %FD
    
    
end
plot(z,utilde,'k--','DisplayName','exact');
legend;
figure;
loglog(2.^[1:Nmax]',error,'r.-');