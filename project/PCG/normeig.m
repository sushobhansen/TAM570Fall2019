function [Sx,Lx]=normeig(Ahx,Bx)

    [Sx,Lx] = eig(full(Ahx),full(Bx));
    
    for i=1:size(Sx,2)
        Sx(:,i) = Sx(:,i)./sqrt(Sx(:,i)'*Bx*Sx(:,i));
    end
    Lx=sparse(Lx);


end