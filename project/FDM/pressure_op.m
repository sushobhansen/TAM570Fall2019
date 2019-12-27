function [output] = pressure_op(direction,pbar,Dhx,Bhx,Bhy,Jpx,Jpy)

if direction==1
    output = Dhx'*Bhx'*(Jpx*pbar*Jpy')*Bhy;
elseif direction==2
    output = Bhx'*(Jpx*pbar*Jpy')*Bhy*Dhx;
else
    output = 'error';
end

end

