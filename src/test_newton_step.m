[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();
delta = 0;
m = 256;
n = 256;

[X1, Y1] = ndgrid(0:m, 1:n); XX1 = X1(:); YY1 = Y1(:);
[X2, Y2] = ndgrid(1:m, 0:n); XX2 = X2(:); YY2 = Y2(:);

f = @(x)exp(-(1 / (1 - ((x / (0.5*m) - 1))^2)));
h1 = arrayfun(f, X1) .* arrayfun(f, Y1);
h2 = arrayfun(f, X2) .* arrayfun(f, Y2);
figure(1), surf(h2)
h1=h1(:);
h2=h2(:);
h = [h1;h2];


ndofb = (m + 1) * n + (n + 1) * m;
ndof = (m - 1) * n + (n - 1) * m;

pl = sparse(1, ndofb);
alpha = sparse(1, ndofb);

P_func = @(pl) P(pl, alpha, XX1, YY1, XX2, YY2, m, n, delta);
t = 0.00001;

d1 = sparse(m + 1, n);d1(2:m-1, 1:n) = zeros(m-2, n);
d1=d1(:);
d2 = sparse(m, n + 1);d2(1:m, 2:m-1) = zeros(m, n-2);
d2(72, 128) = 1.0;
d2=d2(:) / sqrt(d2(:)'*d2(:));

d = [d1;d2];
alpha10_c = alpha(1:(m + 1) * n);
alpha01_c = alpha((m + 1) * n + 1:(m + 1) * n + (n + 1) * m);

val_shift = P_func(h + t*d);
gradPh = (val_shift - P_func(h)) / t;

[gradP_approx, ~, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, h, alpha10_c, alpha01_c, 1.0);

2*gradP_approx' * d - gradPh

function y = P(pl, alpha, XX1, YY1, XX2, YY2, m, n, delta)
ndofb_p1 = (m + 1) * n;
ndofb_p2 = (n + 1) * m;
y = 0;
for i = 1:ndofb_p1
    shift = ndofb_p1;
    tl = i - 1 - (YY1(i) - 1);
    bl = i - (YY1(i) - 1);
    tr = i - 1 - (YY1(i) - 1) + m;
    br = i - (YY1(i) - 1) + m;
    
    if XX1(i) == 0 || XX1(i) == m
        continue
    end
    
    p1 = pl(i);
    p2 = (pl(shift + tr) + pl(shift + br) + pl(shift + tl) + pl(shift + bl)) / 4;
    p_norm = sqrt(p1^2 + p2^2);
    y = y + smooth_pen(p_norm, alpha(i), delta);
end
for i = 1:ndofb_p2
    shift = ndofb_p1;
    tl = i - m - 1 + YY2(i);
    bl = i - m - 1 + YY2(i) + 1;
    tr = i + YY2(i);
    br = i + YY2(i) + 1;
    
    if YY2(i) == 0 || YY2(i) == n
        continue
    end
    
    p1 = pl(i + shift);
    p2 = (pl(tr) + pl(br) + pl(tl) + pl(bl)) / 4;
    
    p_norm = sqrt(p1^2 + p2^2);
    y = y + smooth_pen(p_norm, alpha(i + shift), delta);
end
end

function y = smooth_pen(r, alpha, delta)
    if r >= alpha + delta
        y = 0.5 * (r-alpha)^2 - delta / 2 * (r - alpha) + delta^2 / 6;
    elseif r > alpha && r < alpha + delta 
        y = (r - alpha)^3 / (6 * delta);
    else
        y = 0;
    end
end

