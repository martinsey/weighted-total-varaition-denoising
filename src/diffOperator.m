classdef diffOperator
    %DIFFOPERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        function [d1p, d2p] = fd(n)
            h = 1/n;
            d1p = (eye(n, n - 1) * (-1) + [zeros(1, n - 1); eye(n - 1)]) / h;
            d2p = (eye(n - 1, n) * (-1) + [zeros(n - 1, 1), eye(n - 1)]) / h;
        end
        
        function [d1m, d2m] = bd(n)
            h = 1/n;
            d1m = (eye(n, n - 1) * (-1) + [zeros(1, n - 1); eye(n - 1)]) / h;
            d2m = (eye(n - 1, n) * (-1) + [zeros(n - 1, 1), eye(n - 1)]) / h;
        end
        
        function [d1m, d2m] = bdNeumann(n)
            h = 1/n;
            [d1m, d2m] = diffOperator.bd(n);
            neumann = zeros(n, 1);
            d1m = [neumann, d1m];
            d2m = [neumann.'; d2m];
        end
        
        function [d1p, d2p] = fdNeumann(n)
            h = 1/n;
            [d1p, d2p] = diffOperator.fd(n);
            
            neumann = zeros(n, 1);
            
            d1p = [d1p, neumann];
            d2p = [d2p; neumann.'];
        end
        
        function [ux, uy] = grad(u, n)
            [d1p, d2p] = diffOperator.fd(n);
            ux = u * d1p;
            uy = d2p * u;
        end
        
        function div = div(p1, p2)
            [m, n] = size(p1);
            [d1m, d2m] = diffOperator.bd(n);
            [d1p, d2p] = diffOperator.fd(n);
            div = p1 * d1m + d2m * p2;
        end
    end
end

