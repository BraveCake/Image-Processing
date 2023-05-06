classdef Filters
    methods(Static)
        function result = up_sample(image,rf,cf) %duplicate each column by cf and each row by rf
[r,c,~] = size(image);
result= uint8(zeros(r*rf,c*cf,3));
[i2,~]=deal(1);
for i=1:rf*r
i
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
%method2
function result = down_sample(image,rf,cf) %drop rf-1 rows and cf-1 columns between everytime you move from column to another or row to another
[r,c,~] = size(image);
result= uint8(zeros(floor(r/rf),floor(c/cf),3));
[i2,j2]=deal(1);
for i=0:rf:r-1
    j2=1;
    for j=0:cf:c-1
        result(i2,j2,:)=image(i+1,j+1,:);
        j2=j2+1;
    end
    i2=i2+1;
end
end
%method3
function result = logIntensity(image,c)
result =im2double(image); %converted to double because log expects double valuesubplot(121);
result = im2uint8(c*log(result+1));

end
%method4
function result = reverseLog(image,c) % the reverse function of logIntensity2 it expands dark values and compressed bright ones (makes image darker)
result = uint8(exp(double(image)/c));
end
%method5
function result = powerLawInesity(image,c,gamma) % high gamma = darker image(extends bright) lower gamma = brighter (extends dark)  
result =c*im2double(image).^gamma; 
result = im2uint8(result);

end
%method6
function result = conrastStretching(image,mn,mx) %contrast stretching takes the image, minimum and maximum contrast values 
result= (mx-mn)/(max(image(:))-min(image(:))) * (image-min(image(:))) + mn; % new contrast/old contrast * (pixel-minimum possible value) + new minimum value

end
%method7
function  compareHistogram(img1,img2)
y1 = zeros(256,1);
y2 = zeros(256,1);
[r1,c1,~] = size(img1);
[r2,c2,~] = size(img2);
for(i=1:r1)
    for(j=1:c1)
        y1(img1(i,j)+1)= 1+y1(img1(i,j)+1);
    end
end
for(i=1:r2)
    for(j=1:c2)
        y2(img2(i,j)+1)= 1+y2(img2(i,j)+1);
    end
end
figure("NumberTitle","off","Name","custom");
subplot(1,2,1);
plot([0:255],y1);
title("Original");
subplot(1,2,2);
plot([0:255],y2);
title("Modified");
figure("NumberTitle","off","Name","builtin");
subplot(1,2,1);
imhist(img1);
title("Original");
subplot(1,2,2);
imhist(img2);
title("Modified");
end
%method8
function result= bitSlice(image,b)

result = bitand(2^(b-1),double(image));
end
%method9
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
%method10
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
%method11
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
%method12
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
%method13
function result = medianFilter(image,n)
[r,c] = size(image);
result = image;
temp = zeros(1,n*n);
for(i=1:r)
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
%method14
function result = maxFilter(image,n)
[r,c] = size(image);
result = zeros(r,c);
temp = zeros(n*n,1);
for i=1:r
    i
    for j=1:c
        z = 1;
        for x=max(1,i-floor(n/2)):min(r,i+floor(n/2))
            for y=max(1,j-floor(n/2)):min(c,j+floor(n/2))
                temp(z) = image(x,y);
                z = z + 1;
            end
        end
        result(i,j) = max(temp(:));
    end
end
end
%method15
function result = minFilter(image,n) 
[r,c] = size(image);
result = zeros(r,c);
temp = zeros(n*n,1);
for i=1:r
    i;
    for j=1:c
        z = 1;
        for x=max(1,i-floor(n/2)):min(r,i+floor(n/2))
            for y=max(1,j-floor(n/2)):min(c,j+floor(n/2))
                temp(z) = image(x,y);
                z = z + 1;
            end
        end
        result(i,j) = min(temp(:));
    end
end
end

%method16
function result = laplacian(image,composite) % 2nd derivative of the image if composite is true it uses composite laplace 
% which subtracts the 2nd derivative from original image leading to highlighted edgeds
image=double(image);
if(nargin<2)
    composite=false;
end
[r,c] = size(image);
for(i=2:r-1)
    for(j=2:c-1)
        if(~composite)
        result(i,j)= (image(i,j-1)+image(i,j+1)+image(i+1,j)+image(i-1,j)- 4*image(i,j));
        else
        result(i,j)= 5*double(image(i,j))- image(i,j-1)-image(i,j+1)-image(i+1,j)-image(i-1,j) ;
        end
    end
end
result = uint8(result);

    

end
%method17
function result = robert(image) %first derivative using robert operators (noise sensitive poor detection)
image= double(image);

result = image;
[r,c] = size(image);
for(i=1:r-1)
    for(j=1:c-1)
        result(i,j)= (image(i+1,j+1)-image(i,j)) + ( image(i+1,j)-image(i,j+1)); %applied the following mask
        %-1 0  +  0 -1
        % 0 1     1  0
    end
end
result =uint8(result);
end
%method18
function result = sobel(image) %first derivative using sobel operators (less  noise sensitivity and better detection)
image= double(image);
result = image;
[r,c] = size(image);
for(i=2:r-1)
    for(j=2:c-1)
        x= (-image(i-1,j-1) - 2* image(i-1,j)  - image(i,j+1)) +(image(i+1,j-1) + 2* image(i+1,j)  + image(i+1,j+1)); 
        y = (-image(i-1,j-1) - 2* image(i,j-1)  - image(i+1,j-1)) +(image(i-1,j+1) + 2* image(i,j+1)  + image(i+1,j+1)); 
        result(i,j) = x+y;
    end
end
result =uint8(result);

end

%method(s)19
function compareFourier(img,img2)
FT_img = fft2(img);
FT_img2= fft2(img2);
FT_img = log(abs(FT_img)+1);
FT_img2 = log(abs(FT_img2)+1);
shiftedImg = Filters.myfftshift(FT_img);
shiftedImg2 = Filters.myfftshift(FT_img2);
subplot(1,2,1);
imshow(shiftedImg,[],'initialMagnification','fit')
title('original');
subplot(1,2,2);
imshow(shiftedImg2,[],'initialMagnification','fit')
title('modified');
end
function Fshifted = myfftshift(F)
% Apply an FFT shift to the Fourier transform F of an image
[M,N] = size(F);
p = floor((M+1)/2);
q = floor((N+1)/2);
Fshifted = zeros(M,N);
Fshifted(1:p,1:q) = F(M-p+1:M,N-q+1:N);
Fshifted(1:p,q+1:N) = F(M-p+1:M,1:N-q);
Fshifted(p+1:M,1:q) = F(1:M-p,N-q+1:N);
Fshifted(p+1:M,q+1:N) = F(1:M-p,1:N-q);
end




%method(s)20
function result = ILPF(img,D0)
[M, N] = size(img);
FT_img = fft2(double(img));

% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);
  
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = double(D <= D0);
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
result = uint8(real(ifft2(double(G))));
end
function result = IHPF(img,D0)
[M, N] = size(img);
FT_img = fft2(double(img));

% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);
  
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = double(D >= D0);
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
result = uint8(real(ifft2(double(G))));
end
function result = BLPF(img,D0)
[M, N] = size(img);
FT_img = fft2(double(img));

% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
n= 2*2;
  % Calculating Euclidean Distance
 D = sqrt(U.^2+V.^2); 
 D = D./ D0;
 % Comparing with the cut-off frequency and 
 % determining the filtering mask
 H = 1./((1+D).^n);
%   
% Convolution between the Fourier Transformed
% image and the mask
 G = H.*FT_img;
  
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
result = uint8(real(ifft2(double(G))));
end
function result = BHPF(img,D0)
[M, N] = size(img);
FT_img = fft2(double(img));

% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
n= 2*2;
  % Calculating Euclidean Distance
 D = sqrt(U.^2+V.^2); 
 D = D0./D;
 % Comparing with the cut-off frequency and 
 % determining the filtering mask
 H = 1./((1+D).^n);
%   
% Convolution between the Fourier Transformed
% image and the mask
 G = H.*FT_img;
  
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
result = uint8(real(ifft2(double(G))));
end
function result = GLPF(img, D0)
    [M, N] = size(img);
    FT_img = fft2(double(img));
    
    % Designing filter
    u = 0:(M-1);
    idx = find(u > M/2);
    u(idx) = u(idx) - M;
    v = 0:(N-1);
    idy = find(v > N/2);
    v(idy) = v(idy) - N;
    
    % MATLAB library function meshgrid(v, u) returns
    % 2D grid which contains the coordinates of vectors
    % v and u. Matrix V with each row is a copy 
    % of v, and matrix U with each column is a copy of u
    [V, U] = meshgrid(v, u);
      
    % Calculating Euclidean Distance
    D = sqrt(U.^2 + V.^2);
      
    % Gaussian Lowpass Filter
    H = exp(-(D.^2) / (2 * D0^2));
    
    % Convolution between the Fourier Transformed
    % image and the filter
    G = H .* FT_img;
    
    % Getting the resultant image by Inverse Fourier Transform
    % of the convoluted image using MATLAB library function 
    % ifft2 (2D inverse fast Fourier transform)  
    result = uint8(real(ifft2(double(G))));
end

function result = GHPF(img, D0)
    [M, N] = size(img);
    FT_img = fft2(double(img));
    
    % Designing filter
    u = 0:(M-1);
    idx = find(u > M/2);
    u(idx) = u(idx) - M;
    v = 0:(N-1);
    idy = find(v > N/2);
    v(idy) = v(idy) - N;
    
    % MATLAB library function meshgrid(v, u) returns
    % 2D grid which contains the coordinates of vectors
    % v and u. Matrix V with each row is a copy 
    % of v, and matrix U with each column is a copy of u
    [V, U] = meshgrid(v, u);
      
    % Calculating Euclidean Distance
    D = sqrt(U.^2 + V.^2);
      
    % Gaussian Highpass Filter
    H = 1 - exp(-(D.^2) / (2 * D0^2));
    
    % Convolution between the Fourier Transformed
    % image and the filter
    G = H .* FT_img;
    
    % Getting the resultant image by Inverse Fourier Transform
    % of the convoluted image using MATLAB library function 
    % ifft2 (2D inverse fast Fourier transform)  
    result = uint8(real(ifft2(double(G))));
end
    end %end of static methods
end