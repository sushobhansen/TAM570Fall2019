function output = buoyancy_op(rhoL,rho_0,Bhx,Bhy,g,NE,N)
    output=timestep_semop(rhoL,Bhx,Bhy,NE,N)*(g/rho_0);
end

