originalImage= imread("image.jpg");
figure("NumberTitle","off","Name","Original");
imshow(originalImage);
figure("NumberTitle","off","Name","downsample");
downsampledImage= down_sample(originalImage,4,4);
figure("NumberTitle","off","Name","upsample");
upsampledImage= up_sample(originalImage,3,3);
imshow(upsampledImage);
title("up sampled image");

function result = up_sample(image,rf,cf) %duplicate each column by cf and each row by rf
[r,c,~] = size(image);
result= uint8(zeros(r*rf,c*cf,3));
[i2,~]=deal(1);
i2+"a"
for i=1:rf*r
    j2=1;
    for j=1:cf*c
        result(i,j,:)=image(i2,j2,:);
       if rem(j,cf)==0
           j2=j2+1;
       end
        
    end
    if rem(i,rf)==0
           i2=i2+1;
    end
    
end

end
function result = down_sample(image,rf,cf) %drop rf-1 rows and cf-1 columns between everytime you move from column to another or row to another
[r,c,~] = size(image);
result= uint8(zeros(floor(r/rf),floor(c/cf),3));
[i2,j2]=deal(1);
for i=0:rf:r-1
    j2=1
    for j=0:cf:c-1
        result(i2,j2,:)=image(i+1,j+1,:);
        j2=j2+1;
    end
    i2=i2+1;
end
end
