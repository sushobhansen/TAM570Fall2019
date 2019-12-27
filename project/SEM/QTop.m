function u=QTop(uL,N,NE,Npt,loc_to_glob)
    u=zeros(Npt,1);
    for ie=1:NE
        for i=1:(N+1)*(N+1)
            u(loc_to_glob(i,ie))=u(loc_to_glob(i,ie))+uL(i,ie);
        end
    end
end