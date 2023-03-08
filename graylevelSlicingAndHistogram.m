originalImage = imread('image2.png');
grayImage= rgb2gray(originalImage);
figure("NumberTitle","off","Name","BitSlicing and Histogram")
slicedImage = graySlice(grayImage,8);
for(i=1:8)
slicedImage = graySlice(grayImage,i);
subplot(2,4,i);
imshow(slicedImage);
title("Sliced bit number: "+i);
end

figure("NumberTitle","off","Name","Before equalization Histogram")
subplot(1,2,1)
imshow(grayImage);
subplot(1,2,2)
imhist((grayImage));
figure("NumberTitle","off","Name","After equalization Histogram")
subplot(1,2,1)
outputImage = histogramEqualization(grayImage)
imshow(outputImage);
subplot(1,2,2)
imhist(outputImage);
function result = histogramEqualization(image)
[M,N,~] = size(image);
vector = reshape(image,M*N,1);% 2d to 1d
frequency = histc(vector,[0:255]); %function histc counts the frequency of values within range [1:255]
%- alternative:
%frequency= zeros(256);
%for(i=1:M*N) frequency(vector(i)) = frequency(vector(i))+1; 
% end
%-
sk = round(255/(M*N) * cumsum(frequency)); %L-1/(M*N) * Sgima(Cummlative sum) of frequencies 
result = uint8(changem(image,sk,([0:255]))); %map the pixel values in the image to sk (new values calulated with contrast stretching)
end
function result= graySlice(image,b)

result = bitand(2^(b-1),double(image));

end

