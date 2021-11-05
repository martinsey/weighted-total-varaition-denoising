classdef diffOperator
    %DIFFOPERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
        function [d1p, d2p] = fd(n)
            h = 1/n;
            % grid initialization
            % n = 4;
            %[X, Y] = meshgrid(1:n);
            % omega_0_0 = reshape([X, Y], [n, n, 2]);
    
            % diff-operators
            d1p = (eye(n, n - 1) + [zeros(1, n - 1);(-1) * eye(n - 1)]) / h;
            d2p = (eye(n - 1, n) + [zeros(n - 1, 1), (-1) * eye(n - 1)]) / h;
        end
        
        function [d1m, d2m] = bd(n)
            h = 1/n;
            % grid initialization
            % n = 4;
            %[X, Y] = meshgrid(1:n);
            % omega_0_0 = reshape([X, Y], [n, n, 2]);
    
            % diff-operators
            d1m = (eye(n, n - 1)*(-1) + [zeros(1, n - 1);eye(n - 1)]) / h
            d2m = (eye(n - 1, n)*(-1) + [zeros(n - 1, 1), eye(n - 1)]) / h;
        end
    end
end

