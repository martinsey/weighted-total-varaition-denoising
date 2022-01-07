function [divp, p1, p2] = solve_lower_level(f, alpha, p_0)

[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();

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

dx_f = (f(2:m, 1:n) - f(1:m-1, 1:n)) / hx;
dx_f_vec = sparse(m + 1, n);
dx_f_vec(2:m, 1:n) = dx_f;
dx_f_vec = dx_f_vec(:);

dy_f = (f(1:m, 2:n) - f(1:m, 1:(n-1))) / hy;
dy_f_vec = sparse(m, n + 1);
dy_f_vec(1:m, 2:n) = dy_f;
dy_f_vec = dy_f_vec(:);


[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);

laplace_c = laplace(XX1, YY1, XX2, YY2, m, n);
id = speye(size(laplace_c));
grad_div_ = grad_div(XX1, YY1, XX2, YY2, m, n);

pl = cat(2, p_0{1}, p_0{2});
pl = pl';

pltilde = cat(2, p_0{1}, p_0{2});

g_func = @(p, p_res, eps_) -beta * laplace_c * p - grad_div_ * p + gamma * p - [dx_f_vec; dy_f_vec] + 1/ eps_ * p_res;

counter = 0
while counter <= 100
    [p_res, d_pdelta] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);

    A = -beta * laplace_c - grad_div_ + gamma * id + 1 / eps * d_pdelta;
    p_res = g_func(pl, p_res, eps);
    
    pldelta = A \ -p_res;
    pl = pl + pldelta;
    
    counter = counter + 1
end

p1 = pl(1:(m + 1) * n);
p2 = pl((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));


p1 = reshape(p1, m + 1, n);

dx_p1 = (p1(2:m + 1, 1:n) - p1(1:m, 1:n)) / hx;
p2 = reshape(p2, m, n + 1);
dy_p2 = (p2(1:m, 2:n + 1) - p2(1:m, 1:n)) / hy;

divp = dx_p1 + dy_p2;

end

    
