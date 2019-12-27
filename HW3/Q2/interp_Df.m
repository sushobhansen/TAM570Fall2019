
function [Dum] = interp_Df(m,Df,un)

n  = size(un,1);

uh = Df*rfft(un);
uh = [ uh; zeros(m-n,1) ];

Dum = irfft(uh);

