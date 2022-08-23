function y = d_smooth_pen(r, alpha, delta)
    if r >= delta + alpha
        y = r - alpha - delta /2;
    elseif r > alpha && r < delta + alpha
        y = 0.5 * (r - alpha)^2 / delta;
    else
        y = 0;
    end
end