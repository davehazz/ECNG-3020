
count = 0;
%preallocate for memory

%the vector that holds the images
str = strings(1,100);

%the vector that holds the location of the
location = strings(1,100);

%for loop to store locations
%{

for x = 1:100
    fprintf('hello\n');
    %sets each frame into the respective array
    str(x)= sprintf(" 'C:\\Users\\unity\\source\\repos\\Year two\\clas\\ECNG 3020 Heart Rate\\frame%d.jpg' ",count);
    disp(str(x));
    count= count+1;

end 

%}





%------------PART 2: RGB Channel------------------ 

location = 'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame0.jpg';
image1 = imread(location);
image1R = image1(:,:,1);
image1G = image1(:,:,2);
image1B = image1(:,:,3);
RegRowSt = 500;
RegRowFin = 900;
RegColSt = 500;
RegColFin = 700;

%getting the rgb colour channel
image1GF = double(image1G);
image1RF = double(image1R);
image1BF = double(image1B);


%-------------Part 3: Select a ROI---------------------

%selecting the region of interest
ROI = image1GF(RegRowSt:RegRowFin,RegColSt:RegColFin);
ROI2 = image1RF(RegRowSt:RegRowFin,RegColSt:RegColFin);
ROI3 = image1BF(RegRowSt:RegRowFin,RegColSt:RegColFin);

%{
%displaying the region of interest
figure(1);imshow(uint8(ROI));
figure(2);imshow(uint8(ROI2));
figure(3);imshow(uint8(ROI3));
%}

%-------------Test Part---------------
%{
location = 'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame0.jpg';
image1G = image1d(:,:,2);
RegRowSt = 500;
RegRowFin = 900;
RegColSt = 500;
RegColFin = 700;
image1GF = double(image1G);
ROI = image1GF(RegRowSt:RegRowFin,RegColSt:RegColFin);
%}


%---------------PART 4: Use Filter---------------
%image_test =  'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame9.jpg';
%image1 = imread(RO1);
%image1d = double(image1);


%%% Averaging Filter
%Task 3: Create a 2D averaging filter of size n
%create a variable with the filter size
AvgFilSize = 50;
%Create a 2D empty shell with ones with the size
AvgFil = ones(AvgFilSize,AvgFilSize);
% Make sure the filter sums up to one (normalization)
AvgFil = AvgFil/sum(AvgFil(:));
% apply the filter to the image using convolution
ImgR = convn(double(ROI2),AvgFil,'same');
ImgG = convn(double(ROI),AvgFil,'same');
ImgB = convn(double(ROI3),AvgFil,'same');

% show the image
figure(7);subplot(2,2,3);imshow(uint8(ImgF)); title('ImgF-10');
figure(7);subplot(2,2,4);imshow(uint8(ROI)); title('-original');



%-----------PART 5: Independent Component Analysis---------------

%-----------PART 5, a)signal source creation---------------

%i am finding the mean of the image
% i am doing this to obtain a time series of means of the frames
meanIntensityValueR = mean2(ImgR);
meanIntensityValueG = mean2(ImgG);
meanIntensityValueB = mean2(ImgB);


%{
%comp = 2;
%train ICA model
y = fft2(image1GF);
y = abs(y);
figure(2);
y = log(y+1);

imshow(y,[]);
title('fourier transform of an image')

Mld = rica(y,2);

%plot(Mld.)
%apply ICA
%data_ICA = transform(Mld,)
%}


