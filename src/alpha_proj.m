function [alpha_k] = alpha_proj(alpha_tilde, alpha_up, alpha_down)
%ALPHA_PROJ Execute projection of alpha in H1 onto [alpha_down, alpha_up]
eps_alpha = 1e-12;
tol_p = 1e-4;
size_ = numel(alpha_tilde);
alpha_k = ones(size_, 1)  * alpha_down;
laplace_n_ = laplace_n(sqrt(size_), sqrt(size_));


res = @(alpha) (speye(size_) - laplace_n_)*(alpha - alpha_tilde) + 1 / eps_alpha * (arrayfun(@(v)smooth_max(v, 0), alpha - alpha_up) - arrayfun(@(v)smooth_max(v, 0), alpha_down -  alpha));
h1_dual = @(v) 1 / sqrt(size_) * sqrt(v' * ((speye(size_) - laplace_n_) \ v));

res_0 = res(alpha_k);

while true
    res_k = res(alpha_k);
    xhi_l = arrayfun(@(v)v > alpha_up ||  v  <  alpha_down, alpha_k);
    delta_alpha = (speye(size_) - laplace_n_ + 1 / eps_alpha * spdiags(xhi_l,0, size_, size_)) \ -res_k;
    
    alpha_k = alpha_k + delta_alpha;
    
    fprint("newpton step for executing alpha projection results in error %f", h1_dual(res_k))
    if h1_dual(res_k) < tol_p * h1_dual(res_0)
        break
    end
end

function y = smooth_max(r, delta)
    if r >= delta
        y = r - delta / 2;
    elseif r > 0 && r < delta
        y = 0.5 * r^2 / delta;
    else
        y = 0;
    end
end

end

