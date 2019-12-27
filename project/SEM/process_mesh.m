function [Loc_to_GlobData,CoorPt,NE]=process_mesh(MeshData,N,tol)

tmp0=MeshData(MeshData(:,2)==102,1);
NEpt=MeshData(1,1);
%NEedge=MeshData(1,2);
NEedge=length(tmp0);
NE=MeshData(1,2)-length(tmp0);


CoorEpt=MeshData(2:1+NEpt,2:4);
EEdge=MeshData(2+NEpt:1+NEpt+NEedge,3:4);
Elem=MeshData(2+NEpt+NEedge:end,3:6);
EdgeL=zeros(2,4,NE);
EdgeLgll=zeros(N+1,4,NE);
intPtLgll=zeros(N-1,N-1,NE);
for ie=1:NE
    tmp=[Elem(ie,:),Elem(ie,1)];
    for j=1:4
        EdgeL(:,j,ie)=tmp(j:j+1);
    end
end
EdgeLsort=sort(EdgeL);

CoorPt=CoorEpt;
[z,w]=zwgll(N);

for ie=1:NE
    tmp2=zeros(N-1,4);
    for j=1:4
        J1=(CoorEpt(EdgeL(2,j,ie),:)-CoorEpt(EdgeL(1,j,ie),:))/2;
        J2=(CoorEpt(EdgeL(2,j,ie),:)+CoorEpt(EdgeL(1,j,ie),:))/2;
        tmp=z*J1+J2;
        for i=2:N
            for ipt=1:length(CoorPt(:,1))
                if sum(abs(tmp(i,:)-CoorPt(ipt,:)))<tol
                    tmp2(i-1,j)=ipt;
                    break
                elseif ipt==length(CoorPt(:,1))
                    CoorPt=[CoorPt;tmp(i,:)];
                    tmp2(i-1,j)=ipt+1;
                end
            end
        end
    end
    EdgeLgll(:,:,ie)=[EdgeL(:,:,ie);tmp2];
end



for ie=1:NE
    J1=(CoorPt(flip(EdgeLgll(3:N+1,1,ie)),:)-CoorPt(EdgeLgll(3:N+1,3,ie),:))/2;
    J2=(CoorPt(flip(EdgeLgll(3:N+1,1,ie)),:)+CoorPt(EdgeLgll(3:N+1,3,ie),:))/2;
    for j=1:N-1
        tmp=z*J1(j,:)+J2(j,:);
        for i=2:N
            CoorPt=[CoorPt;tmp(i,:)];
            intPtLgll(i-1,j,ie)=length(CoorPt);
        end
    end
end
Loc_to_GlobData=zeros((N+1)*(N+1),NE);

for ie=1:NE
    Loc_to_GlobData(1,ie)=EdgeLgll(1,2,ie);
    Loc_to_GlobData(N+1,ie)=EdgeLgll(2,2,ie);
    Loc_to_GlobData(2:N,ie)=EdgeLgll(3:end,2,ie);
    tmp=flip(EdgeLgll(3:end,1,ie));
    for j=2:N
        Loc_to_GlobData((j-1)*(N+1)+1,ie)=tmp(j-1);
        Loc_to_GlobData(j*(N+1),ie)=EdgeLgll(j+1,3,ie);
    end
    Loc_to_GlobData(N*(N+1)+1,ie)=EdgeLgll(2,4,ie);
    Loc_to_GlobData(N*(N+1)+2:(N+1)*(N+1)-1,ie)=flip(EdgeLgll(3:end,4,ie));
    Loc_to_GlobData((N+1)*(N+1),ie)=EdgeLgll(1,4,ie);
    for j=2:N
        Loc_to_GlobData((j-1)*(N+1)+2:j*(N+1)-1,ie)=flip(intPtLgll(:,j-1,ie));
    end    
end




end