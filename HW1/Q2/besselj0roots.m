function [roots,iters] = besselj0roots(n)
%Calculates the first n roots of J0
%Calculation is performed using the Secant method
%Returns roots and number of iterations required 
%Tolerance = 1e-14

roots = zeros(n,1);
iters = zeros(n,1);
tol = 1e-14;

xi1 = 2; %Initial guess for first roots

for i=1:n
   value = abs(besselj(0,xi1));
   counter = 0;
   
   while value >= tol
       xi2 = xi1 + 0.1;
       xi1 = xi1 - besselj(0,xi1)*(xi1-xi2)/(besselj(0,xi1)-besselj(0,xi2));
       value = abs(besselj(0,xi1));
       counter = counter + 1;
    end
    
    roots(i) = xi1;
    iters(i) = counter;
    xi1 = xi1 + pi;

end

