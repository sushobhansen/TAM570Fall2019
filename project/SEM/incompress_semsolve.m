function [U1,U2,output] = incompress_semsolve(preconflag,dU1star,dU2star,beta0,Bhx,Bhy,Dhx,Dhy,nu,dt,Jpx,Jpy,NE,N,Np,loc_to_glob,ML1,ML2,Npt)
    dU1star=Maskop(Qop(dU1star,N,NE,loc_to_glob),ML1,NE);
    dU2star=Maskop(Qop(dU2star,N,NE,loc_to_glob),ML2,NE);
    b1=zeros((Np+1)^2,NE); b2=zeros((Np+1)^2,NE);
    for ie=1:NE
        tmpU1=reshape(dU1star(:,ie),[N+1,N+1]);
        tmpU2=reshape(dU2star(:,ie),[N+1,N+1]);
        tmpb1=Jpx'*Bhx(:,:,ie)*(Dhx(:,:,ie)*tmpU1)*Bhy(:,:,ie)'*Jpy+Jpx'*Bhx(:,:,ie)*(tmpU2*Dhy(:,:,ie)')*Bhy(:,:,ie)'*Jpy;
        b1(:,ie)=reshape(tmpb1,[(Np+1)^2,1]);
        %b2(:,ie)=reshape(tmpb2,[(Np+1)^2,1]);
    end
    %b=Maskop(b,ML,NE)
    %b=QTop(Maskop(b1,ML1,NE),Np,NE,Npt,loc_to_glob)+QTop(Maskop(b2,ML2,NE),Np,NE,Npt,loc_to_glob);
    b=QTop(b1,Np,NE,Npt,loc_to_glob);
if(~preconflag)
    [output,xflag] = pcg(@multiplier,b,1e-10,10000);
else
    [output,xflag] = pcg(@multiplier,b,1e-10,10000,@preconditioner);
end

    function precondvec = preconditioner(U)
        %Preconditioner - (beta0*B)^(-1)
        UL=Qop(U,N,NE,loc_to_glob);
        UL=Maskop(UL,ML,NE);
        %%%
        precondvec = precondinv(UL,Bhx,Bhy,NE,N);
        %precondvec = timestep_semop(UL,Bhx.^-1,Bhy.^-1,NE,N);
        precondvec = precondvec/beta0;
        %%%
        precondvec = Maskop(precondvec,ML,NE);
        precondvec = QTop(precondvec,N,NE,Npt,loc_to_glob);
    end

    function Hu = multiplier(U)
       UL=Qop(U,N,NE,loc_to_glob);
       UL=Maskop(UL,ML,NE);
       Hu = beta0* timestep_semop(UL,Bhx,Bhy,NE,N);
       Hu = Hu + nu*dt*timestep_semop(UL,Bhx,Ahy,NE,N);
       Hu = Hu + nu*dt*timestep_semop(UL,Ahx,Bhy,NE,N);
       Hu = Maskop(Hu,ML,NE);
       Hu = QTop(Hu,N,NE,Npt,loc_to_glob);
     end
end

