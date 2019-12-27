function [error] = AD1Dsolver(a,c,L,n)
%TAM 570 AD1D
%Inputs:
%a - diffusion coefficient (scalar)
%c - convection rate (speed) (scalar)
%L - length of domain (scalar)
%n - array of the number of points (vector)

%Returns:
%error - array of errors, one for each value in n (vector)

error = zeros(1,size(n,1));

for i=1:size(n,1)
    h = L/n(i);
    n1 = n(i)-1;
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

    error(i) = max(abs(ub-utilde))/max(abs(utilde));
end
end

