originalImage = imread('image.jpg');
grayImage= rgb2gray(originalImage);
figure("NumberTitle","off","Name","BitSlicing and Histogram")
slicedImage = graySlice(grayImage,8);
for(i=1:8)
slicedImage = graySlice(grayImage,i);
subplot(2,4,i);
imshow(slicedImage);
title("Sliced bit number: "+i);
end
figure("NumberTitle","off","Name","Histogram")
imhist(grayImage);

function result= graySlice(image,b)

result = bitand(2^(b-1),double(image));

end

