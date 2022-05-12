f_orig = readImage("../data/cube.png");

if size(f_orig, 3) == 3
    f_orig = rgb2gray(f_orig);
end

f = readImage("../data/cube.png");

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

alpha = ones(size(f)) * 0.06;

alpha10_c = ext_int_x * alpha(:);
alpha01_c = ext_int_y * alpha(:);

[divp, p, A] = solve_lower_level(f, alpha10_c, alpha01_c, XX1, YY1, XX2, YY2);


c = colormap(lines(20))
rgbImage = ind2rgb(round((divp + f) * 20), c);
imshow(rgbImage);
colorbar();
ax = gca;
% Requires R2020a or later
exportgraphics(ax,'../data/rec_circle_color_coded.png','Resolution',300)
imwrite(f + divp, "../data/rec_cube.png")

ssim(f + divp, f_orig)
psnr(f + divp, f_orig)
