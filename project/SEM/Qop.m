function uL=Qop(u,N,NE,loc_to_glob)
    uL=zeros((N+1)*(N+1),NE);
    for ie=1:NE
        for i=1:(N+1)*(N+1)
            uL(i,ie)=u(loc_to_glob(i,ie));
        end
    end
end