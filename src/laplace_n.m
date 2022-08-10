function laplace_n_ = laplace_n(m, n)
%LAPLACE Summary of this function goes here
%   Detailed explanation goes here


[X, Y] = ndgrid(1:m, 1:n);
XX = X(:);
YY = Y(:);

hx = 1/m;
hy = 1/n;

i1=[];j1=[];v1=[];
for i = 1:m*n
    % corners
    if (XX(i) == 1 && YY(i) == 1)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(1 / (hx^2) + 1 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i + 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i + m;
        v1(5*i - 2) = 1 / (hy^2);
        
        continue
    end
        
    if (XX(i) == 1 && YY(i) == n)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(1 / (hx^2) + 1 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i + 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i - m;
        v1(5*i - 2) = 1 / (hy^2);
        
        continue
    end
        
    if (XX(i) == m && YY(i) == 1)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(1 / (hx^2) + 1 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i - 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i + m;
        v1(5*i - 2) = 1 / (hy^2);
        continue
    end
        
    if (XX(i) == m && YY(i) == n)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(1 / (hx^2) + 1 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i - 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i - m;
        v1(5*i - 2) = 1 / (hy^2);
        continue
    end 
        
        % edges 
    if (XX(i) == 1)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(1 / (hx^2) + 2 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i + 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i - m;
        v1(5*i - 2) = 1 / (hy^2);
        
        i1(5*i - 1) = i;
        j1(5*i - 1) = i + m;
        v1(5*i - 1) = 1 / (hy^2);
        
        continue
    end
        
    if (YY(i) == 1)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(2 / (hx^2) + 1 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i - 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i + 1;
        v1(5*i - 2) = 1 / (hx^2);
        
        i1(5*i - 1) = i;
        j1(5*i - 1) = i + m;
        v1(5*i - 1) = 1 / (hy^2);
        
        continue
    end
        
    if (XX(i) == m)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(1 / (hx^2) + 2 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i - 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i - m;
        v1(5*i - 2) = 1 / (hy^2);
        
        i1(5*i - 1) = i;
        j1(5*i - 1) = i + m;
        v1(5*i - 1) = 1 / (hy^2);
        continue
     end
        
     if (YY(i) == n)
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = -(2 / (hx^2) + 1 / (hy^2));
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i + 1;
        v1(5*i - 3) = 1 / (hx^2);
        
        i1(5*i - 2) = i;
        j1(5*i - 2) = i - 1;
        v1(5*i - 2) = 1 / (hx^2);
        
        i1(5*i - 1) = i;
        j1(5*i - 1) = i - m;
        v1(5*i - 1) = 1 / (hy^2);
        
        continue
     end
     
     i1(5*i - 4) = i;
     j1(5*i - 4) = i;
     v1(5*i - 4) = -(2 / (hx^2) + 2 / (hy^2));
        
     i1(5*i - 3) = i;
     j1(5*i - 3) = i + 1;
     v1(5*i - 3) = 1 / (hx^2);
        
     i1(5*i - 2) = i;
     j1(5*i - 2) = i - m;
     v1(5*i - 2) = 1 / (hy^2);
        
     i1(5*i - 1) = i;
     j1(5*i - 1) = i + m;
     v1(5*i - 1) = 1 / (hy^2);
     
     i1(5*i) = i;
     j1(5*i) = i - 1;
     v1(5*i) = 1 / (hx^2);
end

idx2keep_rows    = sum(abs([i1',j1',v1']),2)>0;
laplace_n_ = [i1', j1', v1'];
laplace_n_ = laplace_n_(idx2keep_rows, :);
laplace_n_ = sparse(laplace_n_(:,1), laplace_n_(:,2), laplace_n_(:,3), m * n, m * n);

end

