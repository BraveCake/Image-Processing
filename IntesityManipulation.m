[originalImage,map]= imread("image2.png");
% originalImage = ind2rgb(originalImage,map);
%originalImage = rgb2gray(im2uint8((originalImage)));
subplot(141)
imshow(originalImage);
outputImage= logIntensity(originalImage,5);
subplot(142)
imshow(outputImage);
subplot(143);
outputImage = powerLawInesity(originalImage,1,.6);
imshow(outputImage)
subplot(144)
imshow(conrastStretching(originalImage,50,250));

function result = logIntensity(image,c)
result =im2double(image); %converted to double because log expects double valuesubplot(121);
result = im2uint8(c*log(result+1));

end
function result = powerLawInesity(image,c,gamma) % high gamma = darker image(extends bright) lower gamma = brighter (extends dark)  
result =c*im2double(image).^gamma; 
result = im2uint8(result);

end
function result = conrastStretching(image,mn,mx) %contrast stretching takes the image, minimum and maximum contrast values 
result= (mx-mn)/(max(image(:))-min(image(:))) * (image-min(image(:))) + mn; % new contrast/old contrast * (pixel-minimum possible value) + new minimum value

end