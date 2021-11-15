function poisson
    n = 3;
    f = ones(n, n);
    f_vec = reshape(f, [n * n, 1]);
    map_vec = 1:n*n;
    map_mat = reshape(map_vec, [n, n]);
    
    bdx = sparse(n*n, n*n);
    for i=1:n
        for j=1:n
            if (i == 1)
                continue
            end
            
            bdx(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  -1;
            bdx(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) = 1;
        end
    end
    
    laplace = sparse(n*n, n*n);
    for i=1:n
        for j=1:n  
            if i == 1 && j == 1
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  2;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i + 1, j))) =  -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j + 1))) =  -1;
                
                continue
            end
            
            if i == n && j == n
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  2;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) =  -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j - 1))) =  -1;
                
                continue
            end
            
            if i == 1 && j == n
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  2;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i + 1, j))) =  -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j - 1))) =  -1;
                
                continue
            end
            
            if i == n && j == 1
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  2;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) =  -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j + 1))) =  -1;
                
                continue
            end
            
            if i == 1
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  3;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i + 1, j))) = -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j - 1))) = -1;
    
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j + 1))) = -1;
                
                continue
            end
            
            if i == n 
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  3;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) = -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j - 1))) = -1;
                laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j + 1))) = -1;
                
                continue
            end
            
            if j == 1
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  3;
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i + 1, j))) = -1;
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) = -1;
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j + 1))) = -1; 
               
               continue
            end
            
            if j == n
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  3;
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i + 1, j))) = -1;
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) = -1;
               laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j - 1))) = -1; 
               
               continue
            end
            
            laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j))) =  4;
            laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i + 1, j))) = -1;
            laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i - 1, j))) = -1;
            laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j + 1))) = -1;
            laplace(map_vec(map_mat(i, j)), map_vec(map_mat(i, j - 1))) = -1;
        end
    end
    
    laplace = laplace / n^2 + eye(n * n, n * n);
    
    b = laplace \ f_vec;
    reshape(b, [n, n]);
end