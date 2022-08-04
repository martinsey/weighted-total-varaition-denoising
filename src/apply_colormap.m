% Define the two colormaps.
k = 255
cmap1 = hot(k)
cmap2 = winter(k) 
% Combine them into one tall colormap.
combinedColorMap = [cmap1; cmap2]
% Pick 15 rows at random.
randomRows = randi(size(combinedColorMap, 1), [k, 1])
% Extract the rows from the combined color map.
randomColors = combinedColorMap(randomRows, :)
load cam_bilevel_rec_staggered.mat
image1 = imread("presentation/par_rec_staggered_alpha_staggered.png");
image2 = imread("presentation/par_rec_regular_alpha_regular.png");
image3 = imread("../data/parrotgray.png");
colormap(randomColors)
rgb = ind2rgb(image1, randomColors);
imwrite(rgb, randomColors, "par_rec_staggered_alpha_staggered_cc.png");
rgb = ind2rgb(image2, randomColors);
imwrite(rgb, randomColors, "par_rec_regular_alpha_regular_cc.png");
rgb = ind2rgb(image3, randomColors);
imwrite(rgb, randomColors, "par_cc.png");
