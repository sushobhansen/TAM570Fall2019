function [U1,U2,dp,iter] = incompress_solve(preconflag,Aflag,beta0,dt,x,y,AAx,AAy,BBx,BBy,dU1star,dU2star,U1b,U2b,Rx1,Ry1,Rx2,Ry2,Dhx,Dhy,Bhx,Bhy,Jpx,Jpy,Jvx,Jvy,Jvxi,Jvyi,LL,UU)

    rhs = (Jpx'*Bhx*Dhx)*(Rx1'*dU1star*Ry1+U1b)*(Bhy'*Jpy);
    rhs = rhs + (Jpx'*Bhx)*(Rx2'*dU2star*Ry2+U2b)*(Dhy'*Bhy'*Jpy);
    rhs = reshape(rhs,size(Jpx,2)*size(Jpy,2),1);
    if(~preconflag)
        [dp,xlfag,relres,iter] = pcg(@multiplier,rhs,1e-6,10000);
    else
        [dp,xlfag,relres,iter] = pcg(@multiplier,rhs,1e-6,10000,@preconditioner);
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
        U = Jvxi'*U*Jvyi;
        %U = Rx1*U*Ry1';
        U = reshape(U,size(Rx1,1)*size(Ry1,1),1);
        U = U(1:end-1);
        precondveci = LL\U;
        precondveci = UU\precondveci;
        precondveci = [precondveci;0];
        precondveci = reshape(precondveci, size(Rx1,1),size(Ry1,1));
        %precondveci = Rx1'*precondveci*Ry1;
        precondveci = Jvxi*precondveci*Jvyi';
        precondveci = reshape(precondveci,size(Jvxi,1)*size(Jvyi,1),1);
        precondveci = precondveci - (sum(precondveci)/size(precondveci,1));
        precondvec = precondveci;
    end
    
    function Ep = multiplier(U)
       U = reshape(U,size(Jpx,2),size(Jpy,2));
       Ep = (AAx*U*BBy') + (BBx*U*AAy'); 
       Ep = reshape(Ep,size(Jpx,2)*size(Jpy,2),1);
    end
end

