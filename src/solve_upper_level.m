function [divp, dx, dy] = solve_upper_level()
f=ones(20,20) * 0.5;
alpha=ones(size(f)) * 0.01;
[m, n] = size(f);
p_0 = {sparse(1, (m + 1) * n), sparse(1, (n + 1) * m)};

noise = 0.07*randn(size(f));

f_noise = f + noise;
[divp, dx, dy] = solve_lower_level(f_noise, alpha, p_0);

f_res = f_noise + divp;

full(noise)
full(divp)


figure, imshow(f)
figure, imshow(f_noise)
figure, imshow(f_res)

end
