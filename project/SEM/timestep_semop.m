function output = timestep_semop(UL,Bhx,Bhy,NE,N)
    output=0*UL;
    for ie=1:NE
        Utmp=reshape(UL(:,ie),[N+1,N+1]);
        res = timestep_op(Utmp,Bhx(:,:,ie),Bhy(:,:,ie));
        res = reshape(res,[(N+1)*(N+1),1]);
        output(:,ie)=res;
    end
end