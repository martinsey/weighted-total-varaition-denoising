function y = readImage(filename)
    y  = double(imread(filename)) / 255;
end