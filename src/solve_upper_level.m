function [f_res] = solve_upper_level()
f = readImage("test_image_64_64.png");
if size(f, 3) == 3
    f = rgb2gray(f);
end
alpha=ones(size(f)) * 0.00025;
[m, n] = size(f);
p_0 = {sparse(1, (m + 1) * n), sparse(1, (n + 1) * m)};

noise = 0.05*randn(size(f));

f_noise = f + noise;
[divp, dx, dy] = solve_lower_level(f_noise, alpha, p_0);

f_res = f_noise + divp;

figure, imshow(f)
figure, imshow(f_noise)
figure, imshow(f_res)

end
