function [f_res] = solve_upper_level(n_w, sigma_up, sigma_down, lambda, alpha_up, alpha_down)
f_noise = readImage("noise.jpeg");
f = readImage("test_image.jpeg");
if size(f_noise, 3) == 3
    f_noise = rgb2gray(f_noise);
end
%noise = randn(size(f_noise)) *  0.05;
%f_noise = f_noise + noise;
% alpha= readmatrix("alpha.txt") * 0.35;

[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();

[m, n] = size(f_noise);
hx = 1/m;
hy = 1/n;
p_0 = {sparse(1, (m + 1) * n), sparse(1, (n + 1) * m)};
tau_0 =  0.1;
c=1e-8;
theta_m=0.25;

%noise = 0.05*randn(size(f));
%f_noise = f + noise;

alpha=ones(size(f_noise)) * 0.001;

[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);

w = ones(n_w, n_w);
R = @(v) conv2(full(v.*v), w, "same");
laplace_n_ = laplace_n(m, n);
H1_norm = @(v) 1 / sqrt(m*n) * sqrt(v' * ((speye(m*n) - laplace_n_) * v));
J = @(alpha, divp) 0.5 * norm(arrayfun(@(v)smooth_max(v, 0),R(divp) - sigma_up^2), "fro") / (n*m) + 0.5 * norm(arrayfun(@(v)smooth_max(v, 0),sigma_down^2 - R(divp)), "fro") / (n*m) + lambda / 2 * H1_norm(alpha)^2;

counter = 0;
while counter <= 6
    [divp, pk, A, alpha10_c, alpha01_c] = solve_lower_level(f_noise, alpha, p_0, XX1, YY1, XX2, YY2);
    b = conv2(arrayfun(@(v)smooth_max(v, 0), R(divp) - sigma_up^2) - arrayfun(@(v)smooth_max(v, 0), sigma_down^2 - R(divp)), w, "same");
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
    
    J_prime = -diag(arrayfun(@(r)dd_smooth_pen(r, delta), pk_norm - [alpha01_c; alpha10_c])) * qk;  % does this depend on alpha_10c and alhpha01c?
    J_prime1 = reshape(J_prime(1:(m + 1) * n), m + 1, n);
    J_prime1 = (J_prime1(1: m, 1:n) + J_prime1(2:m + 1, 1:n)) / 2;
    J_prime2 = reshape(J_prime((m + 1) * n + 1:(m + 1) * n + (n + 1) * m), m, n + 1);
    J_prime2 = (J_prime2(1:m, 1:n) + J_prime2(1:m, 2:n + 1)) / 2;
    J_prime = J_prime1(:) + J_prime2(:) + lambda * (speye(n*m, n*m) - laplace_n_) * alpha(:);
    
    J_grad = lambda * (speye(n*m, n*m) - laplace_n_) \ J_prime;
    
    tau_k = tau_0;
    while true
        alpha_updated  =  alpha_proj(alpha(:)  - tau_k * J_grad , alpha_up, alpha_down);
        [divp_updated, ~, ~, ~, ~] = solve_lower_level(f_noise, alpha_updated, p_0, XX1, YY1, XX2, YY2);
        ssim(f, f_noise + divp_updated)
        psnr(f_noise + divp_updated, f)
        
        if J(alpha_updated, divp_updated) >  J(alpha(:), divp) + c * J_prime' *  (alpha_updated - alpha(:)) &&  tau_k >  10e-5
            tau_k =  theta_m * tau_k
        else
            alpha = alpha_updated;
            figure, imshow(reshape(alpha, m, n) / max(alpha))
            figure, imshow(f_noise + divp)
            J(alpha_updated, divp_updated)
            div_p = divp_updated;
            imwrite(f_noise + divp, "reconstructed_image" + num2str(counter) + ".jpg")
            imwrite(reshape(alpha,  m,  n), "alpha" + num2str(counter) + ".jpg")
            break
        end
    end
    
    if tau_k  <=  10e-6
        break
    end
    
    counter = counter + 1
end

f_res = f_noise + divp;

figure, imshow(f_noise)
figure, imshow(f_res)

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
