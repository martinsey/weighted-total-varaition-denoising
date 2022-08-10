function laplace = laplace(m, n)
%LAPLACE Summary of this function goes here
%   Detailed explanation goes here

[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);
hx = 1/m;
hy = 1/n;

i1 = [];
j1 = [];
v1 = [];

i2 = [];
j2 = [];
v2 = [];


for i = 1:((m + 1) * n)
    if (XX1(i) == 0 || XX1(i) == m) % top and bottom boundary
        continue
    end
    
    if YY1(i) == 1 %left boundary
        i1(5 * i - 4) = i;
        j1(5 * i - 4) = i;
        v1(5 * i - 4) = -(2 / (hx)^2 + 2 / hy^2);
        
        i1(5 * i - 1) = i;
        j1(5 * i - 1) = i + m + 1;
        v1(5 * i - 1) = 1 / hy^2;
        
        i1(5 * i - 2) = i;
        j1(5 * i - 2) = i - 1;
        v1(5 * i - 2) = 1 / hx^2;
        
        i1(5 * i - 3) = i;
        j1(5 * i - 3) = i + 1;
        v1(5 * i - 3) = 1 / hx^2;
        continue
    end
    
    if YY1(i) == n %right boundary
        i1(5 * i - 4) = i;
        j1(5 * i - 4) = i;
        v1(5 * i - 4) = -(2 / hx^2 + 2 / hy^2);
        
        i1(5 * i - 1) = i;
        j1(5 * i - 1) = i - m - 1;
        v1(5 * i - 1) = 1 / hy^2;
        
        i1(5 * i - 2) = i;
        j1(5 * i - 2) = i - 1;
        v1(5 * i - 2) = 1 / hx^2;
        
        i1(5 * i - 3) = i;
        j1(5 * i - 3) = i + 1;
        v1(5 * i - 3) = 1 / hx^2;
        continue
    end
    
    i1(5 * i - 4) = i;
    j1(5 * i - 4) = i;
    v1(5 * i - 4) = -(2 / (hx^2) + 2 / (hy^2));
        
    i1(5 * i - 3) = i;
    j1(5 * i - 3) = i + 1;
    v1(5 * i - 3) = 1 / hx^2;
        
    i1(5 * i - 2) = i;
    j1(5 * i - 2) = i - 1;
    v1(5 * i - 2) = 1 / hx^2;
        
    i1(5 * i - 1) = i;
    j1(5 * i - 1) = i + m + 1;
    v1(5 * i - 1) = 1 / hy^2;
        
    i1(5 * i) = i;
    j1(5 * i) = i - m - 1;
    v1(5 * i) = 1 / hy^2;
end

% symmetrise matrix
for i = 1:((m + 1) * n)
    if XX1(i) == 1
        i1(5 * i - 2) = i;
        j1(5 * i - 2) = i - 1;
        v1(5 * i - 2) = 0;
        
        continue
    end
        
    if XX1(i) == n - 1
        i1(5 * i - 3) = i;
        j1(5 * i - 3) = i + 1;
        v1(5 * i - 3) = 0;
        
        continue
    end
end

for i = 1:((n + 1) * m)
    if YY2(i) == 0 || YY2(i) == n % left and right boundary
        continue
    end
    
    if XX2(i) == 1 % top boundary
        i2(5 * i - 4) = i;
        j2(5 * i - 4) = i;
        v2(5 * i - 4) = -(2 / (hx + 0.5)^2 + 2 / (hy^2));
        
        i2(5 * i - 3) = i;
        j2(5 * i - 3) = i + 1;
        v2(5 * i - 3) = 1 / hx^2;
        
        i2(5 * i - 2) = i;
        j2(5 * i - 2) = i - m;
        v2(5 * i - 2) = 1 / hy^2;
        
        i2(5 * i - 1) = i;
        j2(5 * i - 1) = i + m;
        v2(5 * i - 1) = 1 / hy^2;

        continue
    end
    
    if XX2(i) == m %bottom boundary
        i2(5 * i - 4) = i;
        j2(5 * i - 4) = i;
        v2(5 * i - 4) = -(2 / (hx + 0.5)^2 + 2 / (hy^2));
        
        i2(5 * i - 3) = i;
        j2(5 * i - 3) = i - 1;
        v2(5 * i - 3) = 1 / hx^2;
        
        i2(5 * i - 2) = i;
        j2(5 * i - 2) = i - m;
        v2(5 * i - 2) = 1 / hy^2;
        
        i2(5 * i - 1) = i;
        j2(5 * i - 1) = i + m;
        v2(5 * i - 1) = 1 / hy^2;
        continue
    end
    
    
    i2(5 * i - 4) = i;
    j2(5 * i - 4) = i;
    v2(5 * i - 4) = -(2 / (hx^2) + 2 / (hy^2));
    
    i2(5 * i - 3) = i;
    j2(5 * i - 3) = i + 1;
    v2(5 * i - 3) = 1 / hx^2;
    
    i2(5 * i) = i;
    j2(5 * i) = i - 1;
    v2(5 * i) = 1 / hx^2;
        
    i2(5 * i - 2) = i;
    j2(5 * i - 2) = i - m;
    v2(5 * i - 2) = 1 / hy^2;
        
    i2(5 * i - 1) = i;
    j2(5 * i - 1) = i + m;
    v2(5 * i - 1) = 1 / hy^2;
end
    
for i = 1:((n + 1) * m)
    if YY2(i) == 1
        i2(5 * i - 2) = i;
        j2(5 * i - 2) = i - m;
        v2(5 * i - 2) = 0;
    end
        
    if YY2(i) == n - 1
        i2(5 * i - 1) = i;
        j2(5 * i - 1) = i + m;
        v2(5 * i - 1) = 0;
    end
end

%remove empty cols
idx2keep_rows    = sum(abs([i1',j1',v1']),2)>0 ;
laplace1 = [i1', j1', v1'];
laplace1 = laplace1(idx2keep_rows, :);
laplace1 = sparse(laplace1(:,1), laplace1(:,2), laplace1(:,3), (m + 1) * n, (n + 1) * m);

idx2keep_rows    = sum(abs([i2',j2',v2']),2)>0 ;
laplace2 = [i2', j2', v2'];
laplace2 = laplace2(idx2keep_rows, :);
laplace2 = sparse(laplace2(:,1), laplace2(:,2), laplace2(:,3), m * (n + 1), (m + 1) * n);
laplace = [laplace1, sparse((m + 1) * n,  (n + 1) * m); sparse(m * (n + 1), (m + 1) * n), laplace2];
end

