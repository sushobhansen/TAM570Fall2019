function [bdf,ext] = timecoeffs(k)

bdf = deriv_mat(0:k);
bdf = flip(bdf(end,:));

ext = interp_mat(0,1:k);

end

