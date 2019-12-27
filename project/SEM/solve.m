function output = solve(preconflag,b,beta0,Bhx,Bhy,nu,dt,Ahx,Ahy,NE,N,loc_to_glob,ML,Npt)
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

