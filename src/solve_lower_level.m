max_iter=50;

% calc grad f with zeros boundary as it enforces boundary condition on p
dx_f_vec = sparse(m + 1, n);
dx_f_vec(2:m, 1:n) = reshape(dx * f(:), m - 1, n);
dx_f_vec = dx_f_vec(:);

dy_f_vec = sparse(m, n + 1);
dy_f_vec(1:m, 2:n) = reshape(dy * f(:), m, n - 1);
dy_f_vec = dy_f_vec(:);

%init functions
g_func = @(p, p_pen, eps_l) -beta * laplace_c * p - grad_div_ * p + gamma * p - [dx_f_vec; dy_f_vec] + 1/eps_l * p_pen;
h0div_norm = @(v) sqrt(v' * ((speye((m + 1) * n + m * (n + 1), (m + 1) * n + m * (n + 1)) - grad_div_) \ v) * hx * hy);
g_err = @(p, p_pen, eps_l) h0div_norm(g_func(p, p_pen, eps_l));

%starting values
p_0 = {sparse(1, (m + 1) * n), sparse(1, (n + 1) * m)};
pl = cat(2, p_0{1}, p_0{2});
pl = pl';
pltilde = pl;
eps_l = 0.0001;
[pltilde_pen, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pltilde, alpha10_c, alpha01_c);
pl_pen = [];
p_delta = [];
A_k = [];

%iteration variables
iter_counter = 0;
while true
    if size(A_k, 1) == 0
        [pl_pen, d_pdelta, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);
        A_k = -beta * laplace_c - grad_div_ + gamma * id + 1/eps_l * d_pdelta;
    end
    p_res = g_func(pl, pl_pen, eps_l);
    
    pldelta = A_k \ -p_res;
    pl = pl + pldelta;
    [pl_pen, d_pdelta, ~] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pl, alpha10_c, alpha01_c);
        
    A_k = -beta * laplace_c - grad_div_ + gamma * id + 1/eps_l * d_pdelta;
    
    g_err_new = g_err(pl, pl_pen, eps_l);
    g_err_pltilde = g_err(pltilde, pltilde_pen, eps_l);
    
    fprintf("Executing newton for lower level with error %f and best error %f at eps=%e \n", g_err_new, g_err_pltilde, eps_l);
    if eps_l <= eps
        fprintf("Convergence has been achieved with norm of g %f \n", g_err_new);
        break
    elseif g_err_new <= g_err_pltilde || g_err_new <= 1e-8
        p1 = pl(1:(m + 1) * n);
        p2 = pl((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));
        
        pltilde = pl;
        [pltilde_pen, d_pdelta_tilde, pk_norms] = newton_residual(delta, XX1, YY1, XX2, YY2, m, n, pltilde, alpha10_c, alpha01_c);
        A_k = -beta * laplace_c - grad_div_ + gamma * id + 1/eps_l * d_pdelta_tilde;
        A = A_k;
      
        if g_err(pl, pl_pen, eps_l) < tol_l || iter_counter >= max_iter
            iter_counter = 0;
            eps_l = max(theta_eps * eps_l, eps);
        end
    elseif iter_counter > max_iter
        iter_counter = 0;
        eps_l = max(theta_eps * eps_l, eps);
    else
        iter_counter = iter_counter + 1;
    end
end

p1 = pltilde(1:(m + 1) * n);
p2 = pltilde((m + 1) * n + 1:((m + 1) * n + (n + 1) * m));


p1 = reshape(p1, m + 1, n);

dx_p1 = (p1(2:m + 1, 1:n) - p1(1:m, 1:n)) / hx;
p2 = reshape(p2, m, n + 1);
dy_p2 = (p2(1:m, 2:n + 1) - p2(1:m, 1:n)) / hy;

pk = [p1(:); p2(:)];
divp = dx_p1 + dy_p2;    
