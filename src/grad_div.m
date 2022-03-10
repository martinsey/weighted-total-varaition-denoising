function grad_div = grad_div(XX1, YY1, XX2, YY2, m, n)
grad_div = sparse((m + 1) * n + m * (n + 1), m * (n + 1) + (m + 1) * n);

hx = 1/m;
hy = 1/n;

if isfile("grad_div-" + num2str(m) + "-" +  num2str(n) + ".txt")
    matrix = readmatrix("laplace-" + num2str(m) + "-" +  num2str(n) + ".txt");
    grad_div = sparse(matrix(:,1), matrix(:,2), matrix(:, 3), (m + 1) * n + m * (n + 1), (m + 1) * n + m * (n + 1));
    max(grad_div);
    return
end

for i = 1:((m + 1) * n)
    {"gd 1", i}
    if XX1(i) == 0 || XX1(i) == m % top and bottom boundary
        continue
    end
    
    shift_i = i + (m + 1) * n - YY1(i); %maps to (i, j - 1)
    
    grad_div(i, i + 1) = 1 / hx^2;
    grad_div(i, i) = -2 / hx^2;
    grad_div(i, i - 1) = 1 / hx^2;
    grad_div(i, shift_i + m + 1) = 1 / (hx * hy);
    grad_div(i, shift_i + 1) = -1 / (hx * hy);
    grad_div(i, shift_i + m) = -1 / (hx * hy);
    grad_div(i, shift_i) = 1 / (hx * hy);
end

for i = 1:(m * (n + 1))
    {"gd 2", i}
    shift_i = i + (m + 1) * n;
    % i + YY(2) maps to (i + 1, j - 1)
    
    if YY2(i) == 0 || YY2(i) == n
        continue
    end
    
    grad_div(shift_i, i + YY2(i) - 1 + 2 * (m + 1)) = 1 / (hx * hy);
    grad_div(shift_i, i + YY2(i) - 2 + 2 * (m + 1)) = -1 / (hx * hy);
    grad_div(shift_i, i + YY2(i) - 1 + m + 1) = -1 / (hx * hy);
    grad_div(shift_i, i + YY2(i) - 1) = 1 / (hx * hy);
    grad_div(shift_i, shift_i + m) = 1 / hy^2;
    grad_div(shift_i, shift_i) = -2 / hy^2;
    grad_div(shift_i, shift_i - m) = 1 / hy^2;
end

[i1, i2, v] = find(grad_div);
writematrix([i1, i2, v], "grad_div-" + num2str(m) + "-" +  num2str(n) + ".txt");

end

