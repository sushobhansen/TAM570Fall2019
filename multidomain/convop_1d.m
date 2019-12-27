function cu=convop_1d(Ce,QRt,u);
ul=QRt*u; nl=size(ul,1); N1=size(Ce,1); E=nl/N1;
ul=reshape(ul,N1,E);
cu=Ce*ul; 
cu=QRt'*reshape(cu,N1*E,1);
