function ML=createML(N,NE,boundaryData)
    ML=zeros((N+1)^2,(N+1)^2,NE);
    for ie=1:NE
        ML(:,:,ie)=boundaryData(:,:,ie);
    end
end