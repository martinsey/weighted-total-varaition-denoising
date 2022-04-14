function [grad_div_l, dx, dy, div] = grad_div(m, n)

if isfile("../data/grad_div-" + num2str(m) + "-" +  num2str(n) + ".txt")
    matrix = readmatrix("../data/grad_div-" + num2str(m) + "-" +  num2str(n) + ".txt");
    grad_div_l = sparse(matrix(:,1), matrix(:,2), matrix(:, 3), (m + 1) * n + m * (n + 1), (m + 1) * n + m * (n + 1));
    dx = 0;
    dy = 0;
    div = 0;
    return
end

h = 1 / sqrt(m*n);

div = sparse(m * n,  (m + 1) * n + (n + 1) * m);

[X_div_in, Y_div_in] = ndgrid(1:m, 1:n);
XX_div_in = X_div_in(:);
YY_div_in = Y_div_in(:);

for i = 1: m*n
    {"div x", i}
    if XX_div_in(i) == 1
        div(i, i + YY_div_in(i)) = 1;
        continue
    end
    
    if XX_div_in(i) == m
         div(i, i + YY_div_in(i) - 1) = -1;
        continue
    end
    
    div(i, i + YY_div_in(i) - 1) = -1;
    div(i, i + YY_div_in(i)) = 1;
end

for i = 1: m*n
    {"div y", i}
    if YY_div_in(i) == 1
        div(i, (m + 1) * n + i + m) = 1;
        continue
    end
    
    if YY_div_in(i) == n
        div(i, (m + 1) * n + i) = -1;
        continue
    end
    div(i, (m + 1) * n + i) = -1;
    div(i, (m + 1) * n + i + m) = 1;
end

[X_div, Y_div] = ndgrid(1:m, 1:n + 1);
XX_div = X_div(:);
YY_div = Y_div(:);

dx = sparse((m - 1) * n , m * n);
dy = sparse(m * (n - 1) , m * n);
[X_grad, Y_grad] = ndgrid(1:m - 1, 1:n);
XX_grad = X_grad(:);
YY_grad = Y_grad(:);

for i = 1:(m - 1)*n
    {"dx", i}
    dx(i, i + YY_grad(i) - 1) = -1;
    dx(i, i + YY_grad(i)) = 1;
end

[X_grad, Y_grad] = ndgrid(1:m, 1:n - 1);
XX_grad = X_grad(:);
YY_grad = Y_grad(:);

for i = 1:(n - 1) * m
    {"dy", i}
    dy(i, i) = -1;
    dy(i, i + m) = 1;
end

grad_div_ = [dx; dy] * div;

[X_graddiv, Y_graddiv] = ndgrid(0:m, 1:n);
XX_graddiv = X_graddiv(:);
YY_graddiv = Y_graddiv(:);
grad_div_l = sparse((m + 1) * n + (n + 1) * m, (m + 1) * n + (n + 1) * m);

for i=1:(m + 1) * n
    {"bc x", i}
    if XX_graddiv(i) == 0 || XX_graddiv(i) == m
        continue
    end
    grad_div_l(i, :) = grad_div_(i - 2 * YY_graddiv(i) + 1, :);
end

[X_graddiv, Y_graddiv] = ndgrid(1:m, 0:n);
XX_graddiv = X_graddiv(:);
YY_graddiv = Y_graddiv(:);

j = 0;
for i=1:(n + 1) * m
   {"bc y", i}
   if YY_graddiv(i) == 0 || YY_graddiv(i) == n
       j = j + 1;
       continue
   end
   
   grad_div_l(i + (m + 1) * n, :) = grad_div_(i + (m - 1) * n - j, :);
end

grad_div_l = grad_div_l / h^2;

[i1, i2, v] = find(grad_div_l);
writematrix([i1, i2, v], "../data/grad_div-" + num2str(m) + "-" +  num2str(n) + ".txt");

end

