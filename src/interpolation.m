function [ext_int_x,  ext_int_y] = interpolation(m, n)

if isfile("../data/ext_int_x-" + num2str(m) + "-" +  num2str(n) + ".txt") && isfile("ext_int_y-" + num2str(m) + "-" +  num2str(n) + ".txt")
    matrix1 = readmatrix("../data/ext_int_x-" + num2str(m) + "-" +  num2str(n) + ".txt");
    matrix2 = readmatrix("../data/ext_int_y-" + num2str(m) + "-" +  num2str(n) + ".txt");
    ext_int_x = sparse(matrix1(:,1), matrix1(:,2), matrix1(:, 3), (m + 1) * n, m * n);
    ext_int_y = sparse(matrix2(:,1), matrix2(:,2), matrix2(:, 3), (n + 1) * m, m * n);
    return
end

extend_x = ext_x_bottom(m + 1, n) * ext_x_top(m,n);
extend_y = ext_y_right(m, n + 1) * ext_y_left(m, n);

ext_int_x = int_x(m + 2, n) * extend_x;
ext_int_y = int_y(m, n + 2) * extend_y;

[i1, i2, v] = find(ext_int_x);
writematrix([i1, i2, v], "../data/ext_int_x-" + num2str(m) + "-" +  num2str(n) + ".txt")

[i1, i2, v] = find(ext_int_y);
writematrix([i1, i2, v], "../data/ext_int_y-" + num2str(m) + "-" +  num2str(n) + ".txt")

function interpolate = int_x(m, n)
    [X1, Y1] = ndgrid(1:m - 1, 1:n);
    XX1 = X1(:);
    YY1 = Y1(:);

    interpolate = sparse((m - 1) *  n, m * n);

    for i = 1:(m - 1)*n   
        interpolate(i, i + YY1(i) - 1) = 0.5;
        interpolate(i, i + YY1(i)) = 0.5;
    end
end

function interpolate = int_y(m, n)
    [X1, Y1] = ndgrid(1:m, 1:n - 1);
    XX1 = X1(:);
    YY1 = Y1(:);

    interpolate = sparse(m * (n - 1), m * n);

    for i = 1:m*(n - 1)  
        interpolate(i, i) = 0.5;
        interpolate(i, i + m) = 0.5;
    end
end

function extend_top_x = ext_x_top(m, n)
    [X1, Y1] = ndgrid(0:m, 1:n);
    XX1 = X1(:);
    YY1 = Y1(:);

    extend_top_x = sparse((m + 1) *  n, m * n);

    for i = 1:(m + 1)*n
        if XX1(i) == 0
            extend_top_x(i, i - YY1(i) + 1) = 1;
            continue;
        end
    
        extend_top_x(i, i - YY1(i)) = 1;
    end
end

function extend_bottom_x = ext_x_bottom(m, n)
    extend_bottom_x = sparse((m + 1) *  n, m * n);
    [X2, Y2] = ndgrid(1:m + 1, 1:n);
    XX2 = X2(:);
    YY2 = Y2(:);

    for i = 1:(m + 1)*n
        if XX2(i) == m + 1
            extend_bottom_x(i, i - (YY2(i) - 1) * (m - n) -  1) = 1;
            continue;
        end
    
        extend_bottom_x(i, i - (YY2(i) - 1) * (m - n)) = 1;
    end
end

function extend_left_y  = ext_y_left(m, n)
    [X1, Y1] = ndgrid(1:m, 0:n);
    XX1 = X1(:);
    YY1 = Y1(:);

    extend_left_y = sparse(m * (n + 1), m * n);

    for i = 1:(m + 1)*n
        if YY1(i) == 0
            extend_left_y(i, i) = 1;
            continue;
        end
    
        extend_left_y(i, i - m) = 1;
    end
end

function extend_right_y  = ext_y_right(m, n)
    [X1, Y1] = ndgrid(1:m, 1:n + 1);
    XX1 = X1(:);
    YY1 = Y1(:);

    extend_right_y = sparse(m * (n + 1), m * n);

    for i = 1:(n + 1)*m
        if YY1(i) == n + 1
            extend_right_y(i, i - m) = 1;
            continue;
        end
        
        extend_right_y(i, i) = 1;
    end
end

end

