originalImage= imread("image.jpg");
figure("NumberTitle","off","Name","Negative image filter (GrayScale)")
outputImage= 255-convertRGB2Gray(originalImage);
subplot(1,3,1);
imshow(originalImage);
title("Original")
subplot(132);
imshow(convertRGB2Gray(originalImage),[]);
title("GrayScale")
subplot(133);
imshow(outputImage,[]);
title("Negative filter on GrayScale")

%converting from rgb2gray using average of 3 R G B values
function result = convertRGB2Gray(img)
[r, c,~] = size(img);
result = zeros(r, c); 
for i = 1:r
    for j = 1:c
        result(i, j) = uint8(mean(img(i, j, :)));
    end
end
end
