originalImage = imread('noise.jpg');
grayImage= rgb2gray(originalImage);
subplot(1,3,1)
imshow(grayImage);
title("Original");
subplot(1,3,2);
enhancedImage = averageFilter(grayImage,5);
imshow(enhancedImage);
title("Average");
subplot(1,3,3);
enhancedImage = medianFilter(grayImage,21);
imshow(enhancedImage);
title("Median");
function result = medianFilter(image,n)
[r,c] = size(image);
result = image;
temp = zeros(1,n*n);
for(i=1:r)
    i
    for(j=1:c)
          t=1;
        for(x=i-((n-1)/2):i+((n-1)/2))
            for(y=j-((n-1)/2):j+((n-1)/2))
               
                if(x<1 ||y<1 ||x>r ||y>c) % if there is no pixel in the assumed neighbours replicate the default one
                    temp(t)= image(i,j);
                else
                    temp(t)=image(x,y);
                end
                t=t+1;
            end
            
        end
        
        result(i,j) = median(temp);
    end
    
end

end
function result = medianFilter2(image) %exactly same steps only difference we don't duplicate the pixels that lack to neighbours we just ignore them
%you can't choose neighbours they are 3 only
[r,c] = size(image);
result = image;
for(i=2:r-1)
    for(j=2:c-1)
        result(i,j) = median(reshape([image(i-1:i+1, j-1:j+1)],[1,9]));
    end    
  end
end

function result = averageFilter(image,n)
[r,c] = size(image);
result = zeros(r,c);
for i=1:r
    for j=1:c
        temp = 0;
        count = 0;
        for x=max(1,i-floor(n/2)):min(r,i+floor(n/2))
            for y=max(1,j-floor(n/2)):min(c,j+floor(n/2))
                temp = temp + double(image(x,y));
                count = count + 1;
            end
        end
        result(i,j) = temp / count;
    end
end
result = uint8(result);
end

function result = averageFilter2(image) 
[r,c] = size(image);
result = image;
for(i=2:r-1)
    for(j=2:c-1)
        result(i,j) = mean(reshape([image(i-1:i+1, j-1:j+1)],[1,9]));
    end    
  end
end

function result = maxFilter(image) 
[r,c] = size(image);
result = image;
for(i=2:r-1)
    for(j=2:c-1)
        result(i,j) = max(reshape([image(i-1:i+1, j-1:j+1)],[1,9]));
    end    
  end
end

function result = minFilter(image) 
[r,c] = size(image);
result = image;
for(i=2:r-1)
    for(j=2:c-1)
        result(i,j) = min(reshape([image(i-1:i+1, j-1:j+1)],[1,9]));
    end    
  end
end

