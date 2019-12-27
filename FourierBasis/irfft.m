
%
% Inverse Real fft
%
% Written by Paul Fischer
% Date: Oct 10, 2000
%
%
%  See rfft.m for comments
%
%

function [x] = irfft(r);

n  = size(r,1);
n2 = floor((n-1)/2);
m  = n+1;

x=zeros(n,1);
z=zeros(n,1);

z(1) = n*complex(r(1),0);

for k=1:n2;
    z(1+k) = n*complex(r(2*k),-r(2*k+1))/2;
    z(m-k) = n*complex(r(2*k), r(2*k+1))/2;
end;
if mod(n,2)==0, z(n2+2) = n*complex(r(n),0); end;

x = real( ifft(z) );
