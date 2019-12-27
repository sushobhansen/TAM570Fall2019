function  output=Maskop(UL,ML,NE)
    output=0*UL;
    for ie=1:NE
        output(:,ie)=ML(:,:,ie)*UL(:,ie);
    end
end