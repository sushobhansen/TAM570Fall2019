function [output] = advect_op(direction,U1,U2,rho,Dhx,Dhy,Bmx,Bmy,Jux,Juy)

if direction==0
    output = (Jux'*Bmx*((Jux*(U1)*Juy').*(Jux*Dhx*rho*Juy'))*Bmy'*Juy + Jux'*Bmx*((Jux*U2*Juy').*(Jux*rho*Dhy'*Juy'))*Bmy'*Juy);
elseif direction==1
   output = (Jux'*Bmx*((Jux*(U1)*Juy').*(Jux*Dhx*U1*Juy'))*Bmy'*Juy + Jux'*Bmx*((Jux*U2*Juy').*(Jux*U1*Dhy'*Juy'))*Bmy'*Juy);
elseif direction==2
   output = (Jux'*Bmx*((Jux*(U1)*Juy').*(Jux*Dhx*U2*Juy'))*Bmy'*Juy + Jux'*Bmx*((Jux*U2*Juy').*(Jux*U2*Dhy'*Juy'))*Bmy'*Juy);
else
    output = 'error';
end

end

