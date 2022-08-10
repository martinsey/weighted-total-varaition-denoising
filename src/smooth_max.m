function y = smooth_max(r, delta)
    if r >= delta
        y = r - delta / 2;
    elseif r > 0 && r < delta
        y = 0.5 * r^2 / delta;
    else
        y = 0;
    end
end

