circle_ = zeros(128, 128);
i = 1;
r = 20;

while i <= 128
    j = 1;
    while j <= 128
        if abs(i - 64) < r && abs(j - 64) < r
            circle_(i , j) = 1;
        end
        
        j = j + 1;
    end
    i = i + 1;
end

figure, imshow(circle_)
imwrite(circle_, "../data/circle.png")