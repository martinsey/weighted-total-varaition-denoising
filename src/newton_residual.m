function [p_res, d_pdelta, p_norms] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c)
%NEWTON_RESIDUAL Summary of this function goes here
%   Detailed explanation goes here
 
p_res = zeros((m + 1) * n + (n + 1) * m, 1);
p_norms = zeros((m + 1) * n + (n + 1) * m, 1);

i1 = zeros(1,((m + 1) * n) + (n + 1) * m);
j1 = zeros(1,((m + 1) * n) + (n + 1) * m);
v1 = zeros(1,((m + 1) * n) + (n + 1) * m);
 
i1 = zeros(1,((m + 1) * n) * 5 + (n + 1) * m);
j1 = zeros(1,((m + 1) * n) * 5 + (n + 1) * m);
v1 = zeros(1,((m + 1) * n) * 5 + (n + 1) * m);
 
for i = 1:((m + 1) * n)
    if XX1(i) == 0 || XX1(i) == m % top and bottom boundary
        i1(5*i - 4) = i;
        j1(5*i - 4) = i;
        v1(5*i - 4) = 0;
        
        i1(5*i - 3) = i;
        j1(5*i - 3) = i;
        v1(5*i - 3) = 0;
    
        i1(5*i - 2) = i;
        j1(5*i - 2) = i;
        v1(5*i - 2) = 0;
    
        i1(5*i - 1) = i;
        j1(5*i - 1) = i;
        v1(5*i - 1) = 0;
    
        i1(5*i) = i;
        j1(5*i) = i;
        v1(5*i) = 0;
        continue
    end
    
    shift_i = i + (m + 1) * n - (YY1(i) - 1);
        
    px2 = (pl(shift_i) + pl(shift_i - 1) + pl(shift_i + m - 1) + pl(shift_i + m)) / 4;
    px1 = pl(i);
    
    dpx = dp_val(px1, px2, alpha10_c(i), delta);
    [ddxx, ~] = grad_dp(px1, px2, 0, alpha10_c(i), delta);
 
    i1(5*i - 4) = i;
    j1(5*i - 4) = i;
    v1(5*i - 4) = ddxx;
    
    i1(5*i - 3) = i;
    j1(5*i - 3) = shift_i;
    if (YY2(shift_i - (m + 1) * n) == 0)
        v1(5*i - 3) = 0;
    else
        [~, ddxy] = grad_dp(px1, px2, pl(shift_i), alpha10_c(i), delta);
        v1(5*i - 3) = 1/4 * ddxy;
    end
    
    i1(5*i - 2) = i;
    j1(5*i - 2) = shift_i - 1;
    
    if (YY2(shift_i - 1 - (m + 1) * n) == 0)
        v1(5*i - 2) = 1/4 * 0;
    else
        [~, ddxy] = grad_dp(px1, px2, pl(shift_i - 1), alpha10_c(i), delta);
        v1(5*i - 2) = 1/4 * ddxy;
    end
    
    i1(5*i - 1) = i;
    j1(5*i - 1) = shift_i + m - 1;
    if (YY2(shift_i + m - 1 - (m + 1) * n) == n)
        v1(5*i - 1) = 0;
    else
        [~, ddxy] = grad_dp(px1, px2, pl(shift_i + m - 1), alpha10_c(i), delta);
        v1(5*i - 1) = 1/4 * ddxy;
    end
    
    i1(5*i) = i;
    j1(5*i) = shift_i + m;
    if (YY2(shift_i + m - (m + 1) * n) == n)
        v1(5*i) = 0;
    else
        [~, ddxy] = grad_dp(px1, px2, pl(shift_i + m), alpha10_c(i), delta);
        v1(5*i) = 1/4 * ddxy;
    end
    
    p_res(i) = dpx;
    p_norms(i) = sqrt(px1^2 + px2^2);
end
 
for i = 1:((n + 1) * m)
    shift_i = i + (m + 1) * n;
    if YY2(i) == 0 || YY2(i) == n
        i1(5*(i + (m + 1) * n) - 4) = shift_i;
        j1(5*(i + (m + 1) * n) - 4) = shift_i;
        v1(5*(i + (m + 1) * n) - 4) = 0;
    
        i1(5*(i + (m + 1) * n) - 3) = shift_i;
        j1(5*(i + (m + 1) * n) - 3) = shift_i;
        v1(5*(i + (m + 1) * n) - 3) = 0;
    
        i1(5*(i + (m + 1) * n) - 2) = shift_i;
        j1(5*(i + (m + 1) * n) - 2) = shift_i;
        v1(5*(i + (m + 1) * n) - 2) = 0;
    
        i1(5*(i + (m + 1) * n) - 1) = shift_i;
        j1(5*(i + (m + 1) * n) - 1) = shift_i;
        v1(5*(i + (m + 1) * n) - 1) = 0;
    
        i1(5*(i + (m + 1) * n)) = shift_i;
        j1(5*(i + (m + 1) * n)) = shift_i;
        v1(5*(i + (m + 1) * n)) = 0;
        continue
    end
        
    px2 = pl(shift_i);
    px1 = (pl(i - m - 1 + YY2(i)) + pl(i - m - 1 + YY2(i) + 1) + pl(i + YY2(i)) + pl(i + YY2(i) + 1)) / 4;
    
    dpy = dp_val(px2, px1, alpha01_c(i), delta);
    [ddyy, ~] = grad_dp(px2, px1, 0, alpha01_c(i), delta);
    
    
    i1(5*(i + (m + 1) * n) - 4) = shift_i;
    j1(5*(i + (m + 1) * n) - 4) = shift_i;
    v1(5*(i + (m + 1) * n) - 4) = ddyy;
    
    i1(5*(i + (m + 1) * n) - 3) = shift_i;
    j1(5*(i + (m + 1) * n) - 3) = i - m - 1 + YY2(i);
    if (XX1(i - m - 1 + YY2(i)) == 0)
        v1(5*(i + (m + 1) * n) - 3) = 0;
    else
        [~, ddxy] = grad_dp(px2, px1, pl(i - m - 1 + YY2(i)), alpha01_c(i), delta);
        v1(5*(i + (m + 1) * n) - 3) = 1/4 * ddxy;
    end
    
    i1(5*(i + (m + 1) * n) - 2) = shift_i;
    j1(5*(i + (m + 1) * n) - 2) = i - m - 1 + YY2(i) + 1;
    if XX1(i - m - 1 + YY2(i) + 1) == m
        v1(5*(i + (m + 1) * n) - 2) = 0;
    else
        [~, ddxy] = grad_dp(px2, px1, pl(i - m - 1 + YY2(i) + 1), alpha01_c(i), delta);
        v1(5*(i + (m + 1) * n) - 2) = 1/4 * ddxy;
    end
   
    i1(5*(i + (m + 1) * n) - 1) = shift_i;
    j1(5*(i + (m + 1) * n) - 1) = i + YY2(i);
    if XX1(i + YY2(i)) == 0
        v1(5*(i + (m + 1) * n) - 1) = 0;
    else
        [~, ddxy] = grad_dp(px2, px1, pl(i + YY2(i)), alpha01_c(i), delta);
        v1(5*(i + (m + 1) * n) - 1) = 1/4 * ddxy;
    end
    
    i1(5*(i + (m + 1) * n)) = shift_i;
    j1(5*(i + (m + 1) * n)) = i + YY2(i) + 1;
    if XX1(i + YY2(i) + 1) == m
        v1(5*(i + (m + 1) * n)) = 0;
    else
        [~, ddxy] = grad_dp(px2, px1, pl(i + YY2(i) + 1), alpha01_c(i), delta);
        v1(5*(i + (m + 1) * n)) = 1/4 * ddxy;
    end
 
    p_res(shift_i) = dpy;
    p_norms(shift_i) = sqrt(px1^2 + px2^2);
end
 
d_pdelta = sparse(i1, j1, v1, (m + 1) * n + (n + 1) * m, (m + 1) * n + (n + 1) * m);
 
function y = dp_val(p1, p2, alpha, delta)
    p_norm = sqrt(p1^2 + p2^2);
    if p_norm <= alpha
        y=0;
        return
    end
    y = p1 / p_norm * d_smooth_pen(p_norm, alpha, delta);
end
 
function [ddxx, ddxy] = grad_dp(p1, p2, p2i, alpha, delta)
    p_norm = sqrt(p1^2 + p2^2);
    if p_norm <= alpha
        ddxx = 0;
        ddxy = 0;
        return
    end
    
    ddxx = p2^2 / p_norm^3 * d_smooth_pen(p_norm, alpha, delta) + p1^2 / p_norm^2 * dd_smooth_pen(p_norm, alpha, delta);
    ddxy = -p1*p2 / p_norm^3 * d_smooth_pen(p_norm, alpha, delta) + p1 * p2 / p_norm^2 * dd_smooth_pen(p_norm, alpha, delta);
end
 
function y = d_smooth_pen(r, alpha, delta)
    if r >= delta + alpha
        y = r - alpha - delta /2;
    elseif r > alpha && r < delta + alpha
        y = 0.5 * (r - alpha)^2 / delta;
    else
        y = 0;
    end
end
 
function y = dd_smooth_pen(r, alpha, delta)
    if r >= delta + alpha
        y = 1;
    elseif r > alpha && r < delta + alpha
        y = (r - alpha) / delta;
    else
        y = 0;
    end
end
 
end
