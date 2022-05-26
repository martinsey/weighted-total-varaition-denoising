load cameraman.mat

f_orig = cam;

if size(f_orig, 3) == 3
    f_orig = rgb2gray(f_orig);
end

f = cam_noisy_01;

[beta, gamma, delta, eps, tol_l, theta_eps] = load_variables();

[m, n] = size(f);
hx = 1/m;
hy = 1/n;

[ext_int_x,  ext_int_y] = interpolation(m, n);

[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);
XX1 = X1(:);
YY1 = Y1(:);
XX2 = X2(:);
YY2 = Y2(:);

alpha = ones(size(f)) * 0.0025;

lambda = 1e-9
n_w = 7;
sigma_down = (0.1)^2*(1-(sqrt(2)/n_w));%0.00798;
sigma_up  = (0.1)^2*(1+(sqrt(2)/n_w));%0.01202;
w = ones(n_w, n_w) / n_w^2;
R = @(v) circ_conv(v.^2, w);
laplace_n_ = laplace_n(m, n);
H1_norm = @(v) 1 / sqrt(m*n) * sqrt(v(:)' * ((speye(m*n) - laplace_n_) * v(:)));
J = @(alpha, divp) 0.5 * sum(sum(arrayfun(@(v)smooth_max(v, 0),(R(divp) - sigma_up)).^2)) / (n*m) + 0.5 * sum(sum(arrayfun(@(v)smooth_max(v, 0),(sigma_down - R(divp))).^2)) / (n*m) + 0.5 * lambda * H1_norm(alpha(:))^2;

alpha10_c = ext_int_x * alpha(:);
alpha01_c = ext_int_y * alpha(:);

[divp, p, A, ~] = solve_lower_level(f, alpha10_c, alpha01_c, XX1, YY1, XX2, YY2);
writematrix(divp, "mock_divp.txt");
writematrix(p, "mock_p");
[i1, i2, v] = find(A);
writematrix([i1, i2, v], "mock_A");
imshow(f + divp)
J(alpha, divp)

ssim(f + divp, f_orig)

function y = smooth_max(r, delta)
    if r >= delta
        y = r - delta / 2;
    elseif r > 0 && r < delta
        y = 0.5 * r^2 / delta;
    else
        y = 0;
    end
end
