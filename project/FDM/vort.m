function [output1] = vort(U1,U2,Dhx,Dhy)

output1 = Dhx*U2 - U1*Dhy';

end

