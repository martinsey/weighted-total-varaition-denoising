function divp = solve_lower_level(alpha, beta, gamma, delta, eps, tol_l, theta_eps, f, p_0)

arguments 
    alpha
    beta
    gamma
    delta
    eps
    tol_l
    theta_eps
    f,
    p_0
end

[m, n] = size(f);

alpha01 = (alpha(1:(m-1), 1:n) + alpha(2:m, 1:n)) / 2; % check this
alpha01_c = sparse(m + 1, n);
alpha01_c(2:m, 1:n) = alpha01;

alpha10 = (alpha(1:m, 1:(n-1)) + alpha(1:m, 2:n)) / 2;
alpha10_c = sparse(m, n + 1);
alpha10_c(1:m, 2:n) = alpha10;

alpha01_c = alpha01_c(:);
alpha10_c = alpha10_c(:);

hx = 1/m;
hy = 1/n;

dx_f = (f(2:m, 1:n) - f(1:(m-1), 1:n)) / hx;
dx_f_vec = sparse(m + 1, n);
dx_f_vec(2:m, 1:n) = dx_f;
dx_f_vec = dx_f_vec(:) / hx;

dy_f = (f(1:m, 1:(n-1)) - f(1:m, 2:n)) / hy;
dy_f_vec = sparse(m, n + 1);
dy_f_vec(1:m, 2:n) = dy_f;
dy_f_vec = dy_f_vec(:) / hy;


[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);

laplace_c = laplace(XX1, YY1, XX2, YY2, m, n);
id = speye(size(laplace_c));
grad_div_ = grad_div(XX1, YY1, XX2, YY2, m, n);
    
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

h1inv = eye(size(grad_div_)) - grad_div_;
dual_h1 = @(v)sqrt((v' * h1inv * v) / (hx * hy)); % check this

pl = cat(2, p_0{1}, p_0{2});
pl = pl';
pltilde = cat(2, p_0{1}, p_0{2});
pltilde = pltilde';

d_pdelta = sparse((m + 1) * n + (n + 1) * m, (m + 1) * n + (n + 1) * m);

counter = 0;

while counter <= 100
    p_res = zeros((m + 1) * n + (n + 1) * m, 1);
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
 
    A = -beta * laplace_c - grad_div_ + gamma * id + 1 / eps * d_pdelta;
    p_res = -beta * laplace_c * pl - grad_div_ * pl + gamma * pl - [dx_f_vec; dy_f_vec] + 1/ eps * p_res;
    
    pldelta = A \ -p_res;
    pl = pl + pldelta;
    
    norm(p_res);
    
    counter = counter + 1;
    if counter == 1000
        break
    end
end

p1 = pl(1:(m + 1) * n);
p2 = pl((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));


dx_p1 = reshape(p1, m + 1, n);

dx_p1 = (dx_p1(2:m+1, 1:n) - dx_p1(1:m, 1:n)) / hx;
dy_p2 = reshape(p2, m, n + 1);
dy_p2 = (dy_p2(1:m, 2:n + 1) - dy_p2(1:m, 1:n)) / hy;

divp = dx_p1 + dy_p2;

[m, n] = max(abs(divp(:)));
[x, y] = ind2sub(size(divp),n);

end

    
