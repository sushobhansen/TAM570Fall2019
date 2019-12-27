
function [um] = interp_f2d(m,un)

n  = size(un,1);

uh = rfft(un);
uh = [ uh; zeros(m-n,n) ];

um = irfft(uh);


