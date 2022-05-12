function [v] = circ_conv(u,K)
[n,m]=size(u);

filter_size = size(K);
if m < filter_size(1) || n < filter_size(1)
    ME = MException('circ_conv', 'Image too small');
    throw(ME)
end

[t,]=size(K);
v=zeros(n,m);

ex=(t-1)/2;

u_circ_padded=[u((n-ex+1):n,(m-ex+1):m) u((n-ex+1):n,:)  u((n-ex+1):n,1:ex); 
                u(:,(m-ex+1):m)        u              u(:,1:ex);
                u(1:ex,(m-ex+1):m)     u(1:ex,:)      u(1:ex,1:ex) ];
  
    for i=1:n
        for j=1:m
            
            v(i,j)=sum(sum(K.*u_circ_padded(i:(2*ex+i),j:(2*ex+j))));
            
            
        end
    end
            

            
end

