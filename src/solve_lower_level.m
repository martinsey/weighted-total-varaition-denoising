function [divp, p, A, alpha10_c, alpha01_c] = solve_lower_level(f, alpha, p_0, XX1, YY1, XX2, YY2, noise)

[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();

[m, n] = size(f);

hx = 1/m;
hy = 1/n;

if size(alpha, 2) == 1
    alpha = reshape(alpha, m, n);
end

alpha01 = (alpha(1:(m-1), 1:n) + alpha(2:m, 1:n)) / 2;
alpha01_c = sparse(m + 1, n);
alpha01_c(2:m, 1:n) = alpha01;

alpha10 = (alpha(1:m, 1:(n-1)) + alpha(1:m, 2:n)) / 2;
alpha10_c = sparse(m, n + 1);
alpha10_c(1:m, 2:n) = alpha10;

alpha01_c = alpha01_c(:);
alpha10_c = alpha10_c(:);

dx_f = (f(2:m, 1:n) - f(1:m-1, 1:n)) / hx;
dx_f_vec = sparse(m + 1, n);
dx_f_vec(2:m, 1:n) = dx_f;
dx_f_vec = dx_f_vec(:);

dy_f = (f(1:m, 2:n) - f(1:m, 1:(n-1))) / hy;
dy_f_vec = sparse(m, n + 1);
dy_f_vec(1:m, 2:n) = dy_f;
dy_f_vec = dy_f_vec(:);

laplace_c = laplace(XX1, YY1, XX2, YY2, m, n);
id = speye(size(laplace_c));
grad_div_ = grad_div(XX1, YY1, XX2, YY2, m, n);

pl = cat(2, p_0{1}, p_0{2});
pl = pl';
[pl_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);

g_func = @(p, p_pen, eps_) -beta * laplace_c * p - grad_div_ * p + gamma * p - [dx_f_vec; dy_f_vec] + 1/ eps_ * p_pen;
h0div_norm = @(v) sqrt(v' * (speye((m + 1) * n + m * (n + 1), (m + 1) * n + m * (n + 1))) * v * hx * hy);

pltilde = pl;
[pltilde_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pltilde, alpha10_c, alpha01_c);

eps_l = 1;

h0divptilde = 0;
while eps_l > eps || h0div_norm(g_func(pl, pl_pen, eps_l)) >= tol_l * h0divptilde
    [pl_pen, d_pdelta] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);

    A = -beta * laplace_c - grad_div_ + gamma * id + 1 / eps_l * d_pdelta;
    p_res = g_func(pl, pl_pen, eps_l);

    pldelta = A \ -p_res;
    pl = pl + pldelta;
    [pl_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);
    
    if h0div_norm(g_func(pl, pl_pen, eps_l)) < tol_l * h0div_norm(g_func(pltilde, pltilde_pen, eps_l))
        eps_l = max(theta_eps * eps_l, eps);
        pltilde = pl;
        [pltilde_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);
        h0divptilde = h0div_norm(g_func(pltilde, pltilde_pen, eps_l));
    elseif eps_l <= eps
        break
    end
end

p1 = pl(1:(m + 1) * n);
p2 = pl((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));


p1 = reshape(p1, m + 1, n);

dx_p1 = (p1(2:m + 1, 1:n) - p1(1:m, 1:n)) / hx;
p2 = reshape(p2, m, n + 1);
dy_p2 = (p2(1:m, 2:n + 1) - p2(1:m, 1:n)) / hy;

p = [p1(:), p2(:)];
divp = dx_p1 + dy_p2;

end

    
