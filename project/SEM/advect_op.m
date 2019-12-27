function output= advect_op(C1,C2,U,Dhx,Dhy,Bmx,Bmy,Jux,Juy,NE,N)
        output=0*U;
%output = (Jux'*Bmx*((Jux*(C1)*Juy').*(Jux*Dhx*U*Juy'))*Bmy'*Juy + Jux'*Bmx*((Jux*C2*Juy').*(Jux*U*Dhy'*Juy'))*Bmy'*Juy);
    for ie=1:NE
        tmpC1=reshape(C1(:,ie),[N+1,N+1])*0+0.5;
        tmpC2=reshape(C2(:,ie),[N+1,N+1])*0-0.5;
        tmpU = reshape(U(:,ie),[N+1,N+1]);
        res=Jux'*Bmx(:,:,ie)*((Jux*(tmpC1*Juy')).*(Jux*Dhx(:,:,ie)*tmpU*Juy'))*Bmy(:,:,ie)'*Juy+Jux'*Bmx(:,:,ie)*((Jux*(tmpC2*Juy')).*(Jux*tmpU*Dhy(:,:,ie)'*Juy'))*Bmy(:,:,ie)'*Juy;
        output(:,ie)=reshape(res,[(N+1)^2,1]);
    end


end

