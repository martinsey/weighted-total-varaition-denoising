function laplace = laplace_n(m, n)
%LAPLACE Summary of this function goes here
%   Detailed explanation goes here

laplace = sparse(m * n, m * n);

[X, Y] = ndgrid(1:m, 1:n);
XX = X(:);
YY = Y(:);

hx = 1/m;
hy = 1/n;

if isfile("laplace_n-" + num2str(m) + "-" +  num2str(n) + ".txt")
    matrix = readmatrix("laplace_n-" + num2str(m) + "-" +  num2str(n) + ".txt");
    laplace = sparse(matrix(:,1), matrix(:,2), matrix(:, 3), m*n, m*n);
else
    for i = 1:m*n
        {"laplace_n 1", i}
        
        % corners
        if (XX(i) == 1 && YY(i) == 1)
            laplace(i, i) = -(1 / (hx^2) + 1 / (hy^2));
            laplace(i, i + 1) = 1 / hx^2;
            laplace(i, i + m ) = 1 / hy^2;
            continue
        end
        
        if (XX(i) == 1 && YY(i) == n)
            laplace(i, i) = -(1 / (hx^2) + 1 / (hy^2));
            laplace(i, i + 1) = 1 / hx^2;
            laplace(i, i - m ) = 1 / hy^2;
            continue
        end
        
        if (XX(i) == m && YY(i) == 1)
            laplace(i, i) = -(1 / (hx^2) + 1 / (hy^2));
            laplace(i, i - 1) = 1 / hx^2;
            laplace(i, i + m ) = 1 / hy^2;
            continue
        end
        
        if (XX(i) == m && YY(i) == n)
            laplace(i, i) = -(1 / (hx^2) + 1 / (hy^2));
            laplace(i, i - 1) = 1 / hx^2;
            laplace(i, i - m ) = 1 / hy^2;
            continue
        end 
        
        % edges 
        if (XX(i) == 1)
            laplace(i, i) = -(1 / (hx^2) + 2 / (hy^2));
            laplace(i, i + 1) = 1 / hx^2;
            laplace(i, i + m) = 1 / hy^2;
            laplace(i, i - m) = 1 / hy^2;
            continue
        end
        
        if (YY(i) == 1)
            laplace(i, i) = -(2 / (hx^2) + 1 / (hy^2));
            laplace(i, i + 1) = 1 / hx^2;
            laplace(i, i - 1) = 1 / hx^2;
            laplace(i, i + m) = 1 / hy^2;
            continue
        end
        
        if (XX(i) == m)
            laplace(i, i) = -(1 / (hx^2) + 2 / (hy^2));
            laplace(i, i - 1) = 1 / hx^2;
            laplace(i, i - m) = 1 / hy^2;
            laplace(i, i + m) = 1 / hy^2;
            continue
        end
        
        if (YY(i) == n)
            laplace(i, i) = -(2 / (hx^2) + 1 / (hy^2));
            laplace(i, i + 1) = 1 / hx^2;
            laplace(i, i - 1) = 1 / hx^2;
            laplace(i, i - m) = 1 / hy^2;
            continue
        end
    
        laplace(i, i) = -(2 / (hx^2) + 2 / (hy^2));
        laplace(i, i + 1) = 1 / hx^2;
        laplace(i, i - 1) = 1 / hx^2;
        laplace(i, i + m) = 1 / hy^2;
        laplace(i, i - m) = 1 / hy^2;
    end
    
    [i1, i2, v] = find(laplace);
    
    writematrix([i1, i2, v], "laplace_n-" + num2str(m) + "-" +  num2str(n) + ".txt")
end
end

