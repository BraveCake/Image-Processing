image1 = rgb2gray(imread("ob1.jpg"));
image2 = rgb2gray(imread("ob2.jpg"));
subplot(1,3,1)
imshow(image1)
subplot(132)
imshow(image2)
subplot(133)
imshow(image1-image2);