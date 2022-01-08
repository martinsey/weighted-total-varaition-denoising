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
    
    dp = d_smooth_max(px1 + px2 - alpha01_c(i)) + d_smooth_max(-(px1 + px2 + alpha01_c(i)));
    d_pdelta(i, i) = dp;
    d_pdelta(i, shift_i) = 1/4 * dp;
    d_pdelta(i, shift_i - 1) = 1/4 * dp;
    d_pdelta(i, shift_i + m - 1) = 1/4 * dp;
    d_pdelta(i, shift_i + m) = 1/4 * dp;
    p_res(i) = smooth_max(px1 + px2 - alpha01_c(i)) - smooth_max(-(px1 + px2 + alpha01_c(i)));
end

for i = 1:((n + 1) * m)
    shift_i = i + (m + 1) * n;
    if YY2(i) == 0 || YY2(i) == n
        continue
    end
        
    px2 = pl(shift_i);
    px1 = (pl(i + YY2(i)) + pl(i + YY2(i) + m + 1) + pl(i + YY2(i) + 1) + pl(i + YY2(i) + 1 + m + 1)) / 4;
        
    dp = d_smooth_max(px1 + px2 - alpha10_c(i)) + d_smooth_max(-(px1 + px2 + alpha10_c(i)));
    d_pdelta(shift_i, shift_i) = dp;
    d_pdelta(shift_i, i + YY2(i)) = 1/4 * dp;
    d_pdelta(shift_i, i + YY2(i) + m + 1) = 1/4 * dp;
    d_pdelta(shift_i, i + YY2(i) + 1) = 1/4 * dp;
    d_pdelta(shift_i, i + YY2(i) + 1 + m + 1) = 1/4 * dp;
    p_res(shift_i) = smooth_max(px1 + px2 - alpha10_c(i)) - smooth_max(-(px1 + px2 + alpha10_c(i)));
end

function y = smooth_max(r)
    if r >= delta
        y = r - delta / 2;
    elseif r > 0 && r < delta
        y = 0.5 * r^2 / delta;
    else
        y = 0;
    end
end

function y = d_smooth_max(r)
    if r >= delta
        y = 1;
    elseif r > 0 && r < delta
        y = r / delta;
    else
        y = 0;
    end
end
end

