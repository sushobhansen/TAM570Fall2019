function [U1,U2,dp] = incompress(dU1star,dU2star,dt,beta0,U1b,U2b,SSx,SSy,DDinv,Jpx,Jpy,Bhx,Bhy,Dhx,Dhy,Rx1,Ry1,Rx2,Ry2)
   
    dp = (Jpx'*Bhx*Dhx)*(Rx1'*dU1star*Ry1+U1b)*(Bhy'*Jpy);
    dp = dp + (Jpx'*Bhx)*(Rx2'*dU2star*Ry2+U2b)*(Dhy'*Bhy'*Jpy);
    dp = SSx*((DDinv).*(SSx'*dp*SSy))*SSy';
    dp = -beta0*dp/(dt);
    
    U1 = Rx1*Dhx'*Bhx'*Jpx*dp*Jpy'*Bhy*Ry1';
    U1 = (Rx1*Bhx*Rx1')^(-1)*U1*((Ry1*Bhy*Ry1')^(-1))';
    U1 = dU1star+(dt/beta0)*U1;
    
    U2 = Rx2*Bhx'*Jpx*dp*Jpy'*Bhy*Dhy*Ry2';
    U2 = (Rx2*Bhx*Rx2')^(-1)*U2*((Ry2*Bhy*Ry2')^(-1))';
    U2 = dU2star+(dt/beta0)*U2;
end

