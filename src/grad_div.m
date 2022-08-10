function [grad_div_l, dx, dy, div] = grad_div(m, n)
%GRAD_DIV Get differential opearator representing grad div in a staggered
%grid in col first order.

h = 1 / sqrt(m*n);
div = sparse(m * n,  (m + 1) * n + (n + 1) * m);

i1 = [];
j1 = [];
v1 = [];

[X_div_in, Y_div_in] = ndgrid(1:m, 1:n);
XX_div_in = X_div_in(:);
YY_div_in = Y_div_in(:);

for i = 1: m*n
    if XX_div_in(i) == 1
        i1(4 * i - 1) = i;
        j1(4 * i - 1) = i + YY_div_in(i);
        v1(4 * i - 1) = 1;
        
        continue
    end
    
    if XX_div_in(i) == m
        i1(4 * i) = i;
        j1(4 * i) = i + YY_div_in(i) - 1;
        v1(4 * i) = -1;
        continue
    end
    
    i1(4 * i) = i;
    j1(4 * i) = i + YY_div_in(i) - 1;
    v1(4 * i) = -1;
    
    i1(4 * i - 1) = i;
    j1(4 * i - 1) = i + YY_div_in(i);
    v1(4 * i - 1) = 1;
end

for i = 1: m*n
    if YY_div_in(i) == 1
        i1(4 * i - 2) = i;
        j1(4 * i - 2) = (m + 1) * n + i + m;
        v1(4 * i - 2) = 1;
        continue
    end
    
    if YY_div_in(i) == n
        i1(4 * i - 3) = i;
        j1(4 * i - 3) = (m + 1) * n + i;
        v1(4 * i - 3) = -1;
        continue
    end
    
    i1(4 * i - 3) = i;
    j1(4 * i - 3) = (m + 1) * n + i;
    v1(4 * i - 3) = -1;
    
    i1(4 * i - 2) = i;
    j1(4 * i - 2) = (m + 1) * n + i + m;
    v1(4 * i - 2) = 1;
end

idx2keep_rows    = sum(abs([i1',j1',v1']),2)>0 ;
div = [i1', j1', v1'];
div = div(idx2keep_rows, :);
div = sparse(div(:,1), div(:,2), div(:,3), m * n,  (m + 1) * n + (n + 1) * m);

[X_grad, Y_grad] = ndgrid(1:m - 1, 1:n);
XX_grad = X_grad(:);
YY_grad = Y_grad(:);

i1 = [];
j1 = [];
v1 = [];

for i = 1:(m - 1)*n
    i1(2*i - 1) = i;
    j1(2*i - 1) = i + YY_grad(i) - 1;
    v1(2*i - 1) = -1;
    
    i1(2*i) = i;
    j1(2*i) = i + YY_grad(i);
    v1(2*i) = 1;
end

idx2keep_rows    = sum(abs([i1',j1',v1']),2)>0;
dx = [i1', j1', v1'];
dx = dx(idx2keep_rows, :);
dx = sparse(dx(:,1), dx(:,2), dx(:,3), (m - 1) * n, m * n);

[X_grad, Y_grad] = ndgrid(1:m, 1:n - 1);
XX_grad = X_grad(:);
YY_grad = Y_grad(:);

i1 = [];
j1 = [];
v1 = [];

for i = 1:(n - 1) * m
    i1(2*i - 1) = i;
    j1(2*i - 1) = i;
    v1(2*i - 1) = -1;
    
    i1(2*i) = i;
    j1(2*i) = i + m;
    v1(2*i) = 1;
end

idx2keep_rows    = sum(abs([i1',j1',v1']),2)>0;
dy = [i1', j1', v1'];
dy = dy(idx2keep_rows, :);
dy = sparse(dy(:,1), dy(:,2), dy(:,3), m * (n - 1) , m * n);

grad_div_ = [dx; dy] * div;

[X_graddiv, Y_graddiv] = ndgrid(0:m, 1:n);
XX_graddiv = X_graddiv(:);
YY_graddiv = Y_graddiv(:);
grad_div_l = sparse((m + 1) * n + (n + 1) * m, (m + 1) * n + (n + 1) * m);

i1 = [];
j1 = [];
for i=1:(m + 1) * n
    if XX_graddiv(i) == 0 || XX_graddiv(i) == m
        continue
    end
    
    i1(i) = i - 2 * YY_graddiv(i) + 1;
    j1(i) = i;
end

idx = sum(abs(i1),1)>0;
i1 = i1(idx);

idx = sum(abs(j1),1)>0;
j1 = j1(idx);

grad_div_l(j1, :) = grad_div_(i1, :);

[X_graddiv, Y_graddiv] = ndgrid(1:m, 0:n);
XX_graddiv = X_graddiv(:);
YY_graddiv = Y_graddiv(:);

j = 0;
i1 = [];
j1 = [];
for i=1:(n + 1) * m
   if YY_graddiv(i) == 0 || YY_graddiv(i) == n
       j = j + 1;
       continue
   end
   
   i1(i) = i + (m - 1) * n - j;
   j1(i) = i + (m + 1) * n;
end

idx = sum(abs(i1),1)>0;
i1 = i1(idx);

idx = sum(abs(j1),1)>0;
j1 = j1(idx);

dx = dx/h;
dy = dy/h;

div = div/h;

grad_div_l(j1, :) = grad_div_(i1, :);
grad_div_l = grad_div_l / h^2;

end

