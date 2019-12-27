function [U1,U2,dp] = incompress_solve(preconflag,beta0,dt,AAx,AAy,BBx,BBy,dU1star,dU2star,U1b,U2b,Rx1,Ry1,Rx2,Ry2,Dhx,Dhy,Bhx,Bhy,Jpx,Jpy)

    rhs = (Jpx'*Bhx*Dhx)*(Rx1'*dU1star*Ry1+U1b)*(Bhy'*Jpy);
    rhs = rhs + (Jpx'*Bhx)*(Rx2'*dU2star*Ry2+U2b)*(Dhy'*Bhy'*Jpy);
    rhs = reshape(rhs,size(Jpx,2)*size(Jpy,2),1);
    
    if(~preconflag)
        [dp,xflag] = pcg(@multiplier,rhs,1e-10,10000);
    else
        [dp,xflag] = pcg(@multiplier,rhs,1e-10,10000,@preconditioner);
    end
    
    dp = reshape(dp,size(Jpx,2),size(Jpy,2));
    dp = -beta0*dp/(dt);
    
    U1 = Rx1*Dhx'*Bhx'*Jpx*dp*Jpy'*Bhy*Ry1';
    U1 = (Rx1*Bhx*Rx1')^(-1)*U1*((Ry1*Bhy*Ry1')^(-1))';
    U1 = dU1star+(dt/beta0)*U1;
    
    U2 = Rx2*Bhx'*Jpx*dp*Jpy'*Bhy*Dhy*Ry2';
    U2 = (Rx2*Bhx*Rx2')^(-1)*U2*((Ry2*Bhy*Ry2')^(-1))';
    U2 = dU2star+(dt/beta0)*U2;
    
    function precondvec = preconditioner(U)
        U = reshape(U,size(Jpx,2),size(Jpy,2));
        %pc = ones(size(Jpx,2),1);
        diagAAx = sparse(diag(diag(AAx)));
        diagAAy = sparse(diag(diag(AAy)));
        diagBBx = sparse(diag(diag(BBx)));
        diagBBy = sparse(diag(diag(BBy)));
        precondvec = kron(diagBBy,diagAAx) + kron(diagAAy,diagBBx);
        %precondvec = pc'*precondvec*pc;
        precondvec = precondvec^(-1);
        precondvec = reshape(diag(precondvec),size(Jpx,2),size(Jpy,2));
        precondvec = precondvec*U;
        precondvec = reshape(precondvec,size(Jpx,2)*size(Jpy,2),1);
    end
    
    function Ep = multiplier(U)
       U = reshape(U,size(Jpx,2),size(Jpy,2));
       Ep = (AAx*U*BBy') + (BBx*U*AAy'); 
       Ep = reshape(Ep,size(Jpx,2)*size(Jpy,2),1);
    end
end

