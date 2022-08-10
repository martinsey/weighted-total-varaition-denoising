function [ext_int_x,  ext_int_y] = interpolation(m, n)
% INTERPOLATION Get matrix for moving from Omega00 to OMEGA10 and OMAGE01
% by extending padding the same values to the added axes and linearly
% between two points.

extend_x = ext_x_bottom(m + 1, n) * ext_x_top(m,n);
extend_y = ext_y_right(m, n + 1) * ext_y_left(m, n);

ext_int_x = int_x(m + 2, n) * extend_x;
ext_int_y = int_y(m, n + 2) * extend_y;

function interpolate = int_x(m, n)
    [X1, Y1] = ndgrid(1:m - 1, 1:n);
    XX1 = X1(:);
    YY1 = Y1(:);
    
    i1 = [];j1=[];v1=[];
    for i = 1:(m - 1)*n
        i1(2*i - 1) = i;
        j1(2*i - 1) = i + YY1(i) - 1;
        v1(2*i - 1) = 0.5;
        
        i1(2*i) = i;
        j1(2*i) = i + YY1(i);
        v1(2*i) = 0.5;
    end
    
    interpolate = sparse(i1, j1, v1, (m - 1) *  n, m * n);
end

function interpolate = int_y(m, n)
    [X1, Y1] = ndgrid(1:m, 1:n - 1);
    XX1 = X1(:);
    YY1 = Y1(:);
    
    i1 = [];j1=[];v1=[];
    for i = 1:m*(n - 1)
        i1(2 * i - 1) = i;
        j1(2 * i - 1) = i;
        v1(2 * i - 1) = 0.5;
        
        i1(2*i) = i;
        j1(2*i) = i + m;
        v1(2*i) = 0.5;
    end
    
    interpolate = sparse(i1, j1, v1, m * (n - 1), m * n);
end

function extend_top_x = ext_x_top(m, n)
    [X1, Y1] = ndgrid(0:m, 1:n);
    XX1 = X1(:);
    YY1 = Y1(:);
    
    i1=[];j1=[];v1=[];
    for i = 1:(m + 1)*n
        if XX1(i) == 0
            i1(i) = i;
            j1(i) = i - YY1(i) + 1;
            v1(i) = 1;
            continue;
        end
        
        i1(i) = i;
        j1(i) = i - YY1(i);
        v1(i) = 1;
    end
    
    extend_top_x = sparse(i1, j1, v1, (m + 1) *  n, m * n);
end

function extend_bottom_x = ext_x_bottom(m, n)
    extend_bottom_x = sparse((m + 1) *  n, m * n);
    [X2, Y2] = ndgrid(1:m + 1, 1:n);
    XX2 = X2(:);
    YY2 = Y2(:);

    i1=[];j1=[];v1=[];
    for i = 1:(m + 1)*n
        if XX2(i) == m + 1
            i1(i) = i;
            j1(i) = i - (YY2(i) - 1) * (m - n) -  1;
            v1(i) = 1;
            continue;
        end
        
        i1(i) = i;
        j1(i) = i - (YY2(i) - 1) * (m - n);
        v1(i) = 1;
    end
    
    extend_bottom_x = sparse(i1, j1, v1, (m + 1) *  n, m * n);
end

function extend_left_y  = ext_y_left(m, n)
    [X1, Y1] = ndgrid(1:m, 0:n);
    XX1 = X1(:);
    YY1 = Y1(:);

    i1=[];j1=[];v1=[];
    for i = 1:(m + 1)*n
        if YY1(i) == 0
            i1(i) = i;
            j1(i) = i;
            v1(i) = 1;
            continue;
        end
    
        i1(i) = i;
        j1(i) = i - m;
        v1(i) = 1;
    end
    
    extend_left_y = sparse(i1, j1, v1, m * (n + 1), m * n);
end

function extend_right_y  = ext_y_right(m, n)
    [X1, Y1] = ndgrid(1:m, 1:n + 1);
    XX1 = X1(:);
    YY1 = Y1(:);

    extend_right_y = sparse(m * (n + 1), m * n);
    i1=[];j1=[];v1=[];
    for i = 1:(n + 1)*m
        if YY1(i) == n + 1
            i1(i) = i;
            j1(i) = i - m;
            v1(i) = 1;
            continue;
        end
        
        i1(i) = i;
        j1(i) = i;
        v1(i) = 1;
    end
    
    extend_right_y = sparse(i1, j1, v1, m * (n + 1), m * n);
end

end

