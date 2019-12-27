function [Tr0full,Tpr0full,r,A] = Tsolver(roots,nquads,range,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = range(1); b = range(2);
n = size(roots,1);
num = zeros(n,1);
den = zeros(n,1);
beta = zeros(n,1);
A = zeros(n,1);

%Obtain GLL nodes and weights
[z,w] = zwgll(nquads);
r = a + (b-a)*(z+1.0)*0.5; %rescale [-1,1] to [a,b]

for k=1:n
    num(k) = -1*w'*(besselj(0,roots(k)*r).*r)*(b-a)*0.5;
    den(k) = w'*((besselj(0,roots(k)*r).^2).*r)*(b-a)*0.5;
    beta(k) = num(k)/den(k);
    A(k) = beta(k)/roots(k); 
end

B = -A.*tanh(roots*L);

Tr0full = (B.*besselj(0,roots*r'));
Tpr0full = (beta.*besselj(0,roots*r'));

end

