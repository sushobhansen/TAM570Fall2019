function output = buoyancy_op(rho,rho_0,Bhx,Bhy,g)

    output=(Bhx*rho*Bhy')*(g/rho_0);

end

