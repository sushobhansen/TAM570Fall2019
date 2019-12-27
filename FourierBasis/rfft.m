% Real fft
%
% Written by Paul Fischer
% Date: Oct 10, 2000
%
%
%   Given x(1:n), this returns the vector r(1:n).
%
%   If
%             m_c                              m_s
%      x(i) = sum  a_k cos( 2 pi k (i-1)/n ) + sum  b_k sin( 2 pi k (i-1)/n )
%             k=0                              k=0
%
%   then the coefficients in r() are given by
%
%   r:
%
%     a0   a1     a2     a3     a4     ...     a_m_c
%             b1     b2     b3     b4     ...     b_m_s
%
%   where
%
%           m_s = m_c = n/2 if n is odd
%
%           m_s = m_c-1, m_c = n/2 if n is even
%
%   To invert and recover x(), use  x=irfft(r);
%

function [r] = rfft(x);

n  = size(x,1);
m  = size(x,2);
n2 = floor(n/2);
n1 = n2-1;
if mod(n,2)==1; n1=n2; end;

a = zeros(n,m);
b = zeros(n,m);
z = fft(x);

a0 = z(1,:)/n;
a(1:n-1,:) =  2*real(z(2:end,:))/n;
b(1:n-1,:) = -2*imag(z(2:end,:))/n;

r=zeros(n,m);
r(1,:) = a0;
r(2:2:n-1,:) =  a(1:n1,:);
r(3:2:n  ,:) =  b(1:n1,:);
if mod(n,2)==0, r(n,:) = a(n2,:)/2; end;


