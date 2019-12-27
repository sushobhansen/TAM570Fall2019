function output = precondinv(UL,Bhx,Bhy,NE,N)
    output=0*UL;
    for ie=1:NE
        Utmp=reshape(UL(:,ie),[N+1,N+1]);
        res = timestep_op(Utmp,sparse(Bhx(:,:,ie))^-1,sparse(Bhy(:,:,ie))^-1);
        res = reshape(res,[(N+1)*(N+1),1]);
        output(:,ie)=res;
    end
end