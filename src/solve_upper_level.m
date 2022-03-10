function [f_res] = solve_upper_level(n_w, lambda, alpha_up, alpha_down)

sigma_down = (0.1)^2*(1-(sqrt(2)/n_w));%0.00798;
sigma_up  = (0.1)^2*(1+(sqrt(2)/n_w));%0.01202;
f = readImage("test_image.png");

if size(f, 3) == 3
    f = rgb2gray(f);
end

f_noise = readImage("noise.png");
%f_noise = f - readmatrix("mock_divp.txt");
[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();

[m, n] = size(f_noise);
hx = 1/m;
hy = 1/n;

[ext_int_x,  ext_int_y] = interpolation(m, n);


tau_0 = 0.1;
c=1e-8;
theta_m=0.25;
theta_p=2;

alpha=ones(size(f_noise)) * 0.0005;

[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);


w = ones(n_w, n_w) / n_w^2;
R = @(v) circ_conv(full(v.*v), w);
laplace_n_ = laplace_n(m, n);
H1_norm = @(v) 1 / sqrt(m*n) * sqrt(v' * ((speye(m*n) - laplace_n_) * v));
J = @(alpha, divp) 0.5 * norm(arrayfun(@(v)smooth_max(v, 0),R(divp) - sigma_up), "fro")^2 / (n*m) + 0.5 * norm(arrayfun(@(v)smooth_max(v, 0),sigma_down - R(divp)), "fro")^2 / (n*m) + 0.5 * lambda * H1_norm(alpha)^2;

counter = 0;
divp_updated = [];
tau_k = tau_0;

alpha10_c = ext_int_x * alpha(:);
alpha01_c = ext_int_y * alpha(:);
while counter <= 1000
    if size(divp_updated, 1) == 0
        [divp, pk, A] = solve_lower_level(f_noise, alpha10_c, alpha01_c, XX1, YY1, XX2, YY2);
        imwrite(f_noise + divp, "initial.png")
        
        ssim(f_noise + divp, f)
        psnr(f_noise + divp, f)
        J(alpha(:), divp)
    end
    
    b = circ_conv(arrayfun(@(v)smooth_max(v, 0), R(divp) - sigma_up) - arrayfun(@(v)smooth_max(v, 0), sigma_down - R(divp)), w);
    b = 2*divp(:).* b(:); % check the times two
    b = reshape(b, m, n);
    
    dx_b = sparse(m + 1, n);
    dx_b(2:m, 1:n) = (b(2:m, 1:n) - b(1:m - 1, 1:n)) / hx;
    
    dy_b = sparse(m, n + 1);
    dy_b(1:m, 2:n) = (b(1:m, 2:n) - b(1:m, 1:n - 1)) / hy;
    
    b = [dx_b(:); dy_b(:)];
    
    qk = A \ b;
    
    pk_norm = zeros((m + 1) * n + (n + 1) * m, 1);
    for i = 1:((m + 1) * n)
        if XX1(i) == 0 || XX1(i) == m % top and bottom boundary
            continue
        end
        
        shift_i = i + (m + 1) * n - (YY1(i) - 1);
        
        px2 = (pk(shift_i) + pk(shift_i - 1) + pk(shift_i + m - 1) + pk(shift_i + m)) / 4;
        px1 = pk(i);
        
        pk_norm(i) = sqrt(px1^2 + px2^2);
    end
    
    
    for i = 1:((n + 1) * m)
        shift_i = i + (m + 1) * n;
        if YY2(i) == 0 || YY2(i) == n
            continue
        end
        
        px2 = pk(shift_i);
        px1 = (pk(i + YY2(i)) + pk(i + YY2(i) + m + 1) + pk(i + YY2(i) + 1) + pk(i + YY2(i) + 1 + m + 1)) / 4;
        
        pk_norm(shift_i) = sqrt(px1^2 + px2^2);
    end
    

    J_prime = -1 / eps * diag(ones(size(pk_norm)) .* pk .* arrayfun(@(r)dd_smooth_pen(r, delta), pk_norm - [alpha01_c; alpha10_c])) * qk;  % does this depend on alpha_10c and alhpha01c
    J_prime1 = ext_int_x' * J_prime(1:(m + 1) * n);
    J_prime2 = ext_int_y' * J_prime((m + 1) * n + 1:(m + 1) * n + (n + 1) * m);
    J_prime = J_prime1 + J_prime2 + lambda * (speye(n*m, n*m) - laplace_n_) * alpha(:);
    J_grad = (speye(n*m, n*m) - laplace_n_) \ J_prime;
    
    while true
        alpha_updated = alpha_proj(alpha(:) - tau_k * J_grad, alpha_up, alpha_down);
        alpha10_c = ext_int_x * alpha_updated(:);
        alpha01_c = ext_int_y * alpha_updated(:);
        [divp_updated, pk, A] = solve_lower_level(f_noise, alpha10_c, alpha01_c, XX1, YY1, XX2, YY2);
        
        [J(alpha_updated, divp_updated), J(alpha(:), divp)]
        if J(alpha_updated, divp_updated) >  J(alpha(:), divp) + c * J_prime'*(alpha_updated(:) - alpha(:))
            tau_k =  theta_m * tau_k
        else
            alpha = alpha_updated;
            divp = divp_updated;
            ssim(f_noise + divp, f)
            psnr(f_noise + divp, f)
            tau_k = theta_p * tau_k;
            imwrite(f_noise + divp_updated, "rec_" + num2str(counter) + ".png")
            imwrite(reshape((alpha - min(alpha(:))) / max(alpha(:)) ,  m,  n), "alpha" + num2str(counter) + ".png")
            break
        end
    end
    
    counter = counter + 1
end

f_res = f_noise + divp;

end

% check removed arguments
function y = dd_smooth_pen(r, delta)
    if r >= delta
        y = 1;
    elseif r > 0 && r < delta
        y = r / delta;
    else
        y = 0;
    end
end

function y = smooth_max(r, delta)
    if r >= delta
        y = r - delta / 2;
    elseif r > 0 && r < delta
        y = 0.5 * r^2 / delta;
    else
        y = 0;
    end
end

function y = d_smooth_max(r, delta)
    if r >= delta
        y = 1;
    elseif r > 0 && r < delta
        y = r / delta;
    else
        y = 0;
    end
end
