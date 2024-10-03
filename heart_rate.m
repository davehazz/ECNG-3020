
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





%------------PART 2------------------ 
%{
%location = 'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame0.jpg';
%image1 = imread(location);
image1R = image1d(:,:,1);
image1G = image1d(:,:,2);
image1B = image1d(:,:,3);
RegRowSt = 500;
RegRowFin = 900;
RegColSt = 500;
RegColFin = 700;

%getting the rgb colour channel
image1GF = double(image1G);
image1RF = double(image1R);
image1BF = double(image1B);


%-------------Part 3---------------------

%selecting the region of interest
ROI = image1GF(RegRowSt:RegRowFin,RegColSt:RegColFin);
ROI2 = image1RF(RegRowSt:RegRowFin,RegColSt:RegColFin);
ROI3 = image1BF(RegRowSt:RegRowFin,RegColSt:RegColFin);

%displaying the region of interest
figure(1);imshow(uint8(ROI));
figure(2);imshow(uint8(ROI2));
figure(3);imshow(uint8(ROI3));
%}




