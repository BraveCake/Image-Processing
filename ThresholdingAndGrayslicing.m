originalImage= imread("image2.png");
grayImage = rgb2gray(originalImage);
[method1Result,method2Result] = graySlice(grayImage,120,130);
subplot(1,3,1);
imshow(method1Result);
title("first approach result");
subplot(1,3,2);
imshow(method2Result);
title("second approach result");
subplot(1,3,2);
imshow(applyThresholding(grayImage,150));
title("thresholding result");

%The difference between gray level slicing and thresholding may not clear because they are quite similiar but is easy to grasp from the demonstration given below
% gray level slicing highlights a specific range of gray levels you have a start and end to that range while:-
% thresholding is a special case of contrast stretching it converts pixels above the threshold to white and below the threshold to dark

function [result1,result2]= graySlice(image,start,ennd) %takes 2:3 parametrs representing the start of target pixel value and the end of it
%returns image after gray level slicing with 2 different approaches
%approach 1 set desired pixels to white and undesired ones to black
%approach 2 highlight desired pixels by either black or white and leave the
%undesired pixels values without any change
[r,c] = size(image);
result1 =zeros(r,c);
result2 = zeros(r,c);
%if you uncomment these lines below you can make this function act in a similar manner to thresholding if given two parameters %(image,start aka k)
%if(nargin==2)
%ennd=start;
%end
for(i=1:r)
    for(j=1:c)
        if(image(i,j)>= start && image(i,j)<=ennd)
        result1(i,j)=255;
        result2(i,j)=0;
        else
            result1(i,j)=0;
            result2(i,j)=image(i,j);
        end
    end
end
end

function [result]= applyThresholding(image,k) % if the value >= k set it to white if it is less than k (threshold) set it to black
[r,c] = size(image);
result =zeros(r,c);
for(i=1:r)
    for(j=1:c)
        if(image(i,j)>=k)
        result(i,j)=255;
        else
            result(i,j)=0;
        end
    end
end
end