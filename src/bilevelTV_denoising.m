function [f_res] = bilevelTV_denoising(f, sigma)

% image metadata
[m, n] = size(f);
hx = 1/m;
hy = 1/n;

if m ~= n
    ME = MException('Image must be square, but dimensions are %f x %f', m, n);
    throw(ME)
end

fprintf("Init variables... \n")
% bilevel parameters
n_w=7;      lambda=1e-9;        alpha_up=1e-2;  alpha_down=1e-8;
tau_0=1.0;  c=1e-8;             theta_m=0.25;   theta_p=2;
sigma_down = (sigma)^2*(1-(sqrt(2)/n_w));%0.00798; 
sigma_up  = (sigma)^2*(1+(sqrt(2)/n_w));%0.01202;
fprintf("Init variables... done \n")

w = ones(n_w, n_w) / n_w^2;

% lower level parameters
beta=1e-4;          gamma=1e-4;     delta_initial=0.005; 
eps_initial=5e-7;   tol_l=0.01;     delta_final=0.001;
eps_final=3e-8;     theta_eps=0.6;

fprintf("Init differntial matrices... \n")
%preload matrics
[extend_interpolate_x,  extend_interpolate_y] = interpolation(m, n);
laplace_c = laplace_d(m, n);
id = speye(size(laplace_c));
[grad_div_, dx, dy, div] = grad_div(m, n);
laplace_n_ = laplace_n(m, n);
fprintf("Init differntial matrices... done \n")

%starting values
alpha = ones(m, n) * 0.0025;
tau_k = tau_0;
eps=eps_initial;
delta=delta_initial;

%grid
[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);

%functions
R = @(v) conv2(v.^2, w, "same");
H1_norm = @(v) 1 / sqrt(m*n) * sqrt(v(:)' * ((speye(m*n) - laplace_n_) * v(:)));
J = @(alpha, divp) 0.5 * sum(sum(arrayfun(@(v)smooth_max(v, 0),(R(divp) - sigma_up)).^2)) / (n*m) + 0.5 * sum(sum(arrayfun(@(v)smooth_max(v, 0),(sigma_down - R(divp))).^2)) / (n*m) + 0.5 * lambda * H1_norm(alpha(:))^2;

%iteration variables
divp_updated = [];
alpha10_c = extend_interpolate_x * alpha(:);
alpha01_c = extend_interpolate_y * alpha(:);
prox = [];
is_final = false;
is_first_final = false;

fprintf("Executing gradient descent for reduced objective functional \n")
while true
    if size(divp_updated, 1) == 0 || (is_final && is_first_final)
        solve_lower_level;
    end
    
    b = conv2(arrayfun(@(v)smooth_max(v, 0), R(divp) - sigma_up) - arrayfun(@(v)smooth_max(v, 0), sigma_down - R(divp)), w, "same");
    b = 2*divp(:).* b(:); % check the times two
    b = reshape(b, m, n);
    
    dx_b = sparse(m + 1, n);
    dx_b(2:m, 1:n) = (b(2:m, 1:n) - b(1:m - 1, 1:n)) / hx;
    
    dy_b = sparse(m, n + 1);
    dy_b(1:m, 2:n) = (b(1:m, 2:n) - b(1:m, 1:n - 1)) / hy;
    
    b = [dx_b(:); dy_b(:)];
    
    qk = A \ b;
   
    J_prime = -1 / eps * spdiags(ones(size(pk_norms)) .* pk .* arrayfun(@(r)dd_smooth_pen(r, 0, delta), pk_norms - [alpha01_c; alpha10_c]), 0, (m + 1) * n + (n + 1) * m, (m + 1) * n + (n + 1) * m) * qk;  % does this depend on alpha_10c and alhpha01c
    J_prime1 =  extend_interpolate_x' * J_prime(1:(m + 1) * n);
    J_prime2 = extend_interpolate_y' * J_prime((m + 1) * n + 1:(m + 1) * n + (n + 1) * m);
    J_prime = J_prime1 + J_prime2 + lambda * (speye(n*m, n*m) - laplace_n_) * alpha(:);
    
    J_grad = (speye(n*m, n*m) - laplace_n_) \ J_prime;
    
    J_old = J(reshape(alpha, m, n), reshape(divp, m, n));
    while true
        alpha_updated = alpha_proj(alpha(:) - tau_k * J_grad, alpha_up, alpha_down);
        alpha10_c = extend_interpolate_x * alpha_updated(:);
        alpha01_c = extend_interpolate_y * alpha_updated(:);
        solve_lower_level
        J_new = J(reshape(alpha_updated, m, n), reshape(divp, m, n));

        if J_new >  J_old + c * J_prime'*(alpha_updated(:) - alpha(:))
            tau_k =  theta_m * tau_k
            fprintf("Insufficient decrease of J from %f to %f and set tau=%f \n", J_old, J_new, tau_k)
        else
            prox = [prox, H1_norm(alpha_updated(:) - alpha(:))];
            alpha = alpha_updated;
            tau_k = theta_p * tau_k;
            fprintf("Sufficient decrease of J from %f to %f and set tau=%f \n", J_old, J_new, tau_k)
            u = f + divp;
            save("../data/output/result.mat", "u", "alpha")
            break
        end
    end
    
    prox_measure = prox(end) / prox(1);
    
    if ~is_final && (J_old - J_new) / J_new < 0.03
        tau_k = 1;
        delta = delta_final;
        eps=eps_final;
        is_final=true;
        is_first_final=true;
    end
    
    if prox_measure < 4e-5
        fprintf("Executing gradient descent for reduced objective functional... done with J=%f \n", J_new)
        break
    end
end

f_res = f_noise + divp;

end
