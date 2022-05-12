function [divp, p, A_l] = solve_lower_level(f, alpha10_c, alpha01_c, XX1, YY1, XX2, YY2)
max_iter = 100;


[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();

[m, n] = size(f);

hx = 1/m;
hy = 1/n;

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
[grad_div_, ~, ~, ~] = grad_div(m, n);


p_0 = {sparse(1, (m + 1) * n), sparse(1, (n + 1) * m)};
pl = cat(2, p_0{1}, p_0{2});
pl = pl';

g_func = @(p, p_pen, eps_) -beta * laplace_c * p - grad_div_ * p + gamma * p - [dx_f_vec; dy_f_vec] + 1/ eps_ * p_pen;
h0div_norm = @(v) sqrt(v' * ((speye((m + 1) * n + m * (n + 1), (m + 1) * n + m * (n + 1))) \ v) * hx * hy);

pltilde = pl;
[pltilde_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pltilde, alpha10_c, alpha01_c);

eps_l = 1.0;

pl_pen = [];
p_delta = [];

iter_counter = 0;
while max_iter > iter_counter
    if size(pl_pen, 1) == 0 || size(p_delta, 1) == 0
        [pl_pen, d_pdelta] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);
    end
    
    A = -beta * laplace_c - grad_div_ + gamma * id + 1 / eps_l * d_pdelta;
    p_res = g_func(pl, pl_pen, eps_l);
    
    pldelta = A \ -p_res;
    pl = pl + pldelta;
    [pl_pen, d_pdelta] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);

    A = -beta * laplace_c - grad_div_ + gamma * id + 1 / eps_l * d_pdelta;
    
    [h0div_norm(g_func(pl, pl_pen, eps_l)), tol_l * h0div_norm(g_func(pltilde, pltilde_pen, eps_l))]
    eps_l
    
    if eps_l <= eps && h0div_norm(g_func(pltilde, pltilde_pen, eps_l)) < 1e-3
        break
    elseif h0div_norm(g_func(pl, pl_pen, eps_l)) < tol_l * h0div_norm(g_func(pltilde, pltilde_pen, eps_l))
        iter_counter = 0;
        eps_l = max(theta_eps * eps_l, eps);
        p1 = pl(1:(m + 1) * n);
        p2 = pl((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));


        p1 = reshape(p1, m + 1, n);

        dx_p1 = (p1(2:m + 1, 1:n) - p1(1:m, 1:n)) / hx;
        p2 = reshape(p2, m, n + 1);
        dy_p2 = (p2(1:m, 2:n + 1) - p2(1:m, 1:n)) / hy;

        divp = dx_p1 + dy_p2;
        figure(1), imshow(f + divp);
        drawnow();
        pltilde = pl;
        A_l = A;
        [pltilde_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);
        h0divptilde = h0div_norm(g_func(pltilde, pltilde_pen, eps_l));
    else
        iter_counter = iter_counter + 1;
    end
    
    iter_counter = iter_counter + 1;
end

p1 = pltilde(1:(m + 1) * n);
p2 = pltilde((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));


p1 = reshape(p1, m + 1, n);

dx_p1 = (p1(2:m + 1, 1:n) - p1(1:m, 1:n)) / hx;
p2 = reshape(p2, m, n + 1);
dy_p2 = (p2(1:m, 2:n + 1) - p2(1:m, 1:n)) / hy;

p = [p1(:); p2(:)];
divp = dx_p1 + dy_p2;

end

    
