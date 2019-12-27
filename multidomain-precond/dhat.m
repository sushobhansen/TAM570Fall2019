      function[D] =  dhatn(x)
%
%     Compute Dhat
%
      
      ni = length(x);
      a  = ones(ni,1);
      for i=1:ni;
        for j=1:(i-1);  a(i)=a(i)*(x(i)-x(j)); end;
        for j=(i+1):ni; a(i)=a(i)*(x(i)-x(j)); end;
      end;
      a=1./a; % These are the alpha_i's

      D=zeros(ni,ni);
      for j=1:ni; for i=1:ni; D(i,j)=x(i)-x(j); end; D(j,j)=1; end;
      D=1./D;
      for i=1:ni; D(i,i)=0; D(i,i)=sum(D(i,:)); end;

      for j=1:ni; for i=1:ni;
         if i~=j; D(i,j) = a(j)/( a(i)*(x(i)-x(j)));end;
      end;end;

