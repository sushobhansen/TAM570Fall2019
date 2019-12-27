function [um] = interp_f(M,un)

[n,m] = size(un);

uh = rfft(un);
uh = [uh;zeros(M-n,1)];
um = irfft(uh);

end

