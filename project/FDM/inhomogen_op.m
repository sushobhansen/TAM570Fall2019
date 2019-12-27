function [output1] = inhomogen_op(beta0,nu,dt,U1b,Ahx,Ahy,Bhx,Bhy)

output1 = (beta0*Bhx*U1b*Bhy' + nu*dt*(Bhx*U1b*Ahy' + Ahx*U1b*Bhy'));

end

