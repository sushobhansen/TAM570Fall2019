format compact; format shorte; close all

N=8;      %% Polynomial Order
M=40+3*N; %% Target. Increase M>>N for smoother output.

%[z,w]=zwglc(N);  % Chebyshev
%[z,w]=zwuni(N);  % Uniform
[z,w]=zwgll(N);  % GLL

[zm,wm] = zwuni(M); J=interp_mat(zm,z); % Fine points on GLL
[R,S] = ndgrid(z,z);                    % Tensor of points on Omega-hat
F=0.0*R; i=2; j=3;  F(i,j) = 1.0;           % Plot phi_i,j

FF=J*F*J';X=J*R*J';Y=J*S*J'; 
hold off; surfl(X,Y,FF,[1 2 3],[.35 .10 .1 01]);   %%% Lighted grayscale surface
colormap('gray'); shading interp; axis on; hold on;
xlabel('r'); ylabel('s'); zlabel('F');
title(sprintf('l_i l_j for i=%d, j=%d',i,j));

%Plot lines
LX=J*F;  RX=J*R;  RY=J*S;  plot3(RX ,RY ,LX ,'k-','linewidth',1.3)
LX=F*J'; RX=R*J'; RY=S*J'; plot3(RX',RY',LX','k-','linewidth',1.3)
axis equal
print -dpng 'cardinal_2d.png'

% Plot derivatives
Dh = deriv_mat(z);
FFp = J*(Dh*F)*J';
figure;
hold off; surfl(X,Y,FFp,[1 2 3],[.35 .10 .1 01]);   %%% Lighted grayscale surface
colormap('gray'); shading interp; axis on; hold on;
xlabel('r'); ylabel('s'); zlabel('F');
title(sprintf('(l_i)_r l_j for i=%d, j=%d',i,j));

% Generate deformed geometry on coarse 2x2 grid
 
 figure;
 for k=0:2
   theta=(2*k*pi)/3 + 2*pi*(0:3)/3; x=cos(theta); y=sin(theta);
   xmid=zeros(3,1);ymid=xmid;
   xmid(:)=0.5*(x(1:3)+x(2:end)); Xc=[ 0 xmid(3); xmid(1) x(1) ]'; 
   ymid(:)=0.5*(y(1:3)+y(2:end)); Yc=[ 0 ymid(3); ymid(1) y(1) ]';
   Nc=1; [zc,wc]=zwgll(Nc);
   Jc=interp_mat(z,zc); X=Jc*Xc*Jc'; Y=Jc*Yc*Jc';
 
   F=0*R;
   if k==0;   F(i,j) = 1; end;      % Plot phi_i,j
   if i*j==1; F(i,j) = 1; end;      % Plot phi_i,j
   FF=J*F*J';XF=J*X*J';YF=J*Y*J'; surf(XF,YF,FF,'EdgeColor','none'); hold on
   LX=J*F;  RX=J*X;  RY=J*Y;  plot3(RX ,RY ,LX ,'k-','linewidth',1.3)
   LX=F*J'; RX=X*J'; RY=Y*J'; plot3(RX',RY',LX','k-','linewidth',1.3)
 end
xlabel('x'); ylabel('y'); zlabel('F'); 
title(sprintf('l_i l_j in Transformed Domain for i=%d, j=%d',i,j));
print -dpng 'midside_2d.png'



