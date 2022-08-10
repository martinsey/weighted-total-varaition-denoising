function y = dd_smooth_pen(r, alpha, delta)
    if r >= delta + alpha
        y = 1;
    elseif r > alpha && r < delta + alpha
        y = (r - alpha) / delta;
    else
        y = 0;
    end
end
