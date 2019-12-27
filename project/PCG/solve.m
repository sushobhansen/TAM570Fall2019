function [output,iter] = solve(preconflag,b,beta0,Bhx,Bhy,nu,dt,Ahx,Ahy,Rx,Ry)
b = reshape(b,size(Rx,1)*size(Ry,1),1);
if(~preconflag)
    [output,xflag,relres,iter] = pcg(@multiplier,b,1e-10,10000);
else
    [output,xflag,relres,iter] = pcg(@multiplier,b,1e-10,10000,@preconditioner);
end

output = reshape(output,size(Rx,1),size(Ry,1));

    function precondvec = preconditioner(U)
        %Preconditioner - (beta0*B)^(-1)
        U = reshape(U,size(Rx,1),size(Ry,1));
        precondvec = (Rx*Bhx*Rx')^(-1)*U*(Ry*Bhy'*Ry')^(-1);
        precondvec = precondvec/beta0;
        precondvec = reshape(precondvec,size(Rx,1)*size(Ry,1),1);
    end

    function Hu = multiplier(U)
       U = reshape(U,size(Rx,1),size(Ry,1));
       Hu = beta0*((Rx*Bhx*Rx')*U*(Ry*Bhy'*Ry'));
       Hu = Hu + nu*dt*((Rx*Ahx*Rx')*U*(Ry*Bhy'*Ry'));
       Hu = Hu + nu*dt*((Rx*Bhx*Rx')*U*(Ry*Ahy'*Ry'));
       
       Hu = reshape(Hu,size(Rx,1)*size(Ry,1),1);
    end
end

