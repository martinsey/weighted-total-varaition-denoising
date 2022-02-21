function laplace = laplace(XX1, YY1, XX2, YY2, m, n)
%LAPLACE Summary of this function goes here
%   Detailed explanation goes here

laplace1 = sparse((m + 1) * n, (m + 1) * n);
laplace2 = sparse(m * (n + 1), m * (n + 1));

hx = 1/m;
hy = 1/n;

if isfile("laplace-" + num2str(m) + "-" +  num2str(n) + ".txt")
    matrix = readmatrix("laplace-" + num2str(m) + "-" +  num2str(n) + ".txt");
    laplace = sparse(matrix(:,1), matrix(:,2), matrix(:, 3), (m + 1) * n + m * (n + 1), (m + 1) * n + m * (n + 1));
else
    for i = 1:((m + 1) * n)
        {"laplace 1", i}
        [XX1(i), YY1(i)]
        if (XX1(i) == 0 || XX1(i) == m) % top and bottom boundary
            continue
        end
    
        if YY1(i) == 1 %left boundary
            laplace1(i, i) = -(2 / hx^2 + 1 / hy^2);
            laplace1(i, i + m + 1) = 1 / hy^2;
            laplace1(i, i - 1) = 1 / hx^2;
            laplace1(i, i + 1) = 1 / hx^2;
            continue
        end
    
        if YY1(i) == n %right boundary
            laplace1(i, i) = -(2 / hx^2 + 1 / hy^2);
            laplace1(i, i - m - 1) = 1 / hy^2;
            laplace1(i, i + 1) = 1 / hx^2;
            laplace1(i, i - 1) = 1 / hx^2;
            continue
        end
        
        if XX1(i) == 1
            laplace1(i, i) = -(2 / (hx^2) + 2 / (hy^2));
            laplace1(i, i + 1) = 1 / hx^2;
            laplace1(i, i + m + 1) = 1 / hy^2;
            laplace1(i, i - m - 1) = 1 / hy^2;
            continue
        end
        
        
        if XX1(i) == m - 1
            laplace1(i, i) = -(2 / (hx^2) + 2 / (hy^2));
            laplace1(i, i - 1) = 1 / hx^2;
            laplace1(i, i + m + 1) = 1 / hy^2;
            laplace1(i, i - m - 1) = 1 / hy^2;
            continue
        end
    
        laplace1(i, i) = -(2 / (hx^2) + 2 / (hy^2));
        laplace1(i, i + 1) = 1 / hx^2;
        laplace1(i, i - 1) = 1 / hx^2;
        laplace1(i, i + m + 1) = 1 / hy^2;
        laplace1(i, i - m - 1) = 1 / hy^2;
    end
    
    for i = 1:((m + 1) * n)
        if XX1(i) == 1
            laplace1(i, i - 1) = 0;
            continue
        end
        
        if XX1(i) == n - 1
            laplace1(i, i + 1) = 0;
            continue
        end
    end

    for i = 1:((n + 1) * m)
        {"laplace 2", i}
        if YY2(i) == 0 || YY2(i) == n % left and right boundary
            continue
        end
    
        if XX2(i) == 1 % top boundary
            laplace2(i, i) = -(1 / hx^2 + 2 / hy^2);
            laplace2(i, i + 1) = 1 / hx^2;
            laplace2(i, i - m) = 1 / hy^2;
            laplace2(i, i + m) = 1 / hy^2;
            continue
        end
    
        if XX2(i) == m %bottom boundary
            laplace2(i, i) = -(1 / hx^2 + 2 / hy^2);
            laplace2(i, i - 1) = 1 / hx^2;
            laplace2(i, i - m) = 1 / hy^2;
            laplace2(i, i + m) = 1 / hy^2;
            continue
        end
    
        laplace2(i, i) = -(2 / (hx^2) + 2 / (hy^2));
        laplace2(i, i + 1) = 1 / hx^2;
        laplace2(i, i - 1) = 1 / hx^2;
        laplace2(i, i + m) = 1 / hy^2;
        laplace2(i, i - m) = 1 / hy^2;
    end
    
    for i = 1:((n + 1) * m)
        if YY2(i) == 1
            laplace2(i, i - m) = 0;
        end
        
        if  YY2(i) == n - 1
            laplace2(i, i + m) = 0;
        end
    end

    laplace = [laplace1, sparse((m + 1) * n,  (n + 1) * m); sparse(m * (n + 1), (m + 1) * n), laplace2];
    
    [i1, i2, v] = find(laplace);
    
    writematrix([i1, i2, v], "laplace-" + num2str(m) + "-" +  num2str(n) + ".txt")
end
end

