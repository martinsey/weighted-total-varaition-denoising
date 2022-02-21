function [p_res, d_pdelta] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c)
%NEWTON_RESIDUAL Summary of this function goes here
%   Detailed explanation goes here

d_pdelta = sparse((m + 1) * n + (n + 1) * m, (m + 1) * n + (n + 1) * m);
p_res = sparse((m + 1) * n + (n + 1) * m, 1);

for i = 1:((m + 1) * n)
    if XX1(i) == 0 || XX1(i) == m % top and bottom boundary
        continue
    end
    
    shift_i = i + (m + 1) * n - (YY1(i) - 1);
        
    px2 = (pl(shift_i) + pl(shift_i - 1) + pl(shift_i + m - 1) + pl(shift_i + m)) / 4;
    px1 = pl(i);
    
    dpx = dpx_val(px1, px2, alpha01_c(i), delta);
    [ddxx, ddxy] = grad1_dp(px1, px2, alpha01_c(i), delta);
    d_pdelta(i, i) = ddxx;
    d_pdelta(i, shift_i) = 1/4 * ddxy;
    d_pdelta(i, shift_i - 1) = 1/4 * ddxy;
    d_pdelta(i, shift_i + m - 1) = 1/4 * ddxy;
    d_pdelta(i, shift_i + m) = 1/4 * ddxy;
    p_res(i) = dpx;
end

for i = 1:((n + 1) * m)
    shift_i = i + (m + 1) * n;
    if YY2(i) == 0 || YY2(i) == n
        continue
    end
        
    px2 = pl(shift_i);
    px1 = (pl(i + YY2(i)) + pl(i + YY2(i) + m + 1) + pl(i + YY2(i) + 1) + pl(i + YY2(i) + 1 + m + 1)) / 4;
    
    dpy = dpy_val(px1, px2, alpha10_c(i), delta);
    [ddyy, ddxy] = grad2_dp(px1, px2, alpha10_c(i), delta);
    d_pdelta(shift_i, shift_i) = ddyy;
    d_pdelta(shift_i, i + YY2(i)) = 1/4 * ddxy;
    d_pdelta(shift_i, i + YY2(i) + m + 1) = 1/4 * ddxy;
    d_pdelta(shift_i, i + YY2(i) + 1) = 1/4 * ddxy;
    d_pdelta(shift_i, i + YY2(i) + 1 + m + 1) = 1/4 * ddxy;
    p_res(shift_i) = dpy;
end

function y = dpx_val(p1, p2, alpha, delta)
    p_norm = sqrt(p1^2 + p2^2);
    if p_norm <= alpha
        y=0;
        return
    end
    y = p1 / norm([p1, p2]) * d_smooth_pen(p_norm, alpha, delta);
end

function y = dpy_val(p1, p2, alpha, delta)
    p_norm = sqrt(p1^2 + p2^2);
    if p_norm <= alpha
        y=0;
        return
    end
    
    y = p2 / norm([p1, p2]) * d_smooth_pen(p_norm, alpha, delta);
end

function [ddxx, ddxy] = grad1_dp(p1, p2, alpha, delta)
    p_norm = sqrt(p1^2 + p2^2);
    
    if p_norm <= alpha
        ddxx = 0;
        ddxy = 0;
        return
    end
    
    ddxx = p2^2 / p_norm^3 * d_smooth_pen(p_norm, alpha, delta) + p1^2 / p_norm^2 * dd_smooth_pen(p_norm, alpha, delta);
    ddxy = p1 * p2 / p_norm^3 * d_smooth_pen(p_norm, alpha, delta) + p1 * p2 / p_norm^2 * dd_smooth_pen(p_norm, alpha, delta);
end

function [ddyy, ddxy] = grad2_dp(p1, p2, alpha, delta)
    p_norm = sqrt(p1^2 + p2^2);
    
    if p_norm <= alpha
        ddyy = 0;
        ddxy = 0;
        return
    end
    
    ddyy = p1^2 / p_norm^3 * d_smooth_pen(p_norm, alpha, delta) + p2^2 / p_norm^2 * dd_smooth_pen(p_norm, alpha, delta);
    ddxy = p1 * p2 / p_norm^3 * d_smooth_pen(p_norm, alpha, delta) + p1 * p2 / p_norm^2 * dd_smooth_pen(p_norm, alpha, delta);
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

