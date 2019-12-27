
%
%  Test block rfft and irfft scripts.
%

%  If u  = u [Nx,Ny], then uh= rfft(u )=uh[Nx,Ny].

%  If uh = uh[Nx,Ny], then u =irfft(uh)=u [Nx,Ny].

%  The routines rfft and irfft allow you to convert data from physical to wave
%  space in the "x" direction (i.e., in the leading dimension of the input data.



Nx=90; Ny=20;
[zs,ws]=zwgll(Ny);
[zr,wr]=zwuni(Nx); zr=pi*(zr+1);
[X,Y]=ndgrid(zr,zs);
u=(1-Y.*Y).*cos(10*X).*sin(X);
uh=rfft(u);
w =irfft(uh);
mesh(X,Y,uh); pause
mesh(X,Y,w)
