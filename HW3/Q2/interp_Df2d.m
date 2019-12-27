
function [Dum] = interp_Df2d(m,Df,un)

n  = size(un,1);

uh = Df*rfft(un);
uh = [ uh; zeros(m-n,n) ];

Dum = irfft(uh);

