function [R_filtered_avg,G_filtered_avg,B_filtered_avg] = filter_img_avg(inputArg1,inputArg2,inputArg3)

%this function performs the filering on the image
R_filtered_avg = zeros(151,301,800);
G_filtered_avg = zeros(151,301,800);
B_filtered_avg = zeros(151,301,800);


%%% Averaging Filter
%Task 3: Create a 2D averaging filter of size n
%create a variable with the filter size
AvgFilSize = 7;
%Create a 2D empty shell with ones with the size
AvgFil = ones(AvgFilSize,AvgFilSize);
% Make sure the filter sums up to one (normalization)
AvgFil = AvgFil/sum(AvgFil(:));
% apply the filter to the image using convolution
%ImgR = convn(double(ROI2),AvgFil,'same');
%ImgG = convn(double(ROI),AvgFil,'same');
%ImgB = convn(double(ROI3),AvgFil,'same');

for x = 1:700
    R_filtered_avg(:,:,x) = convn(inputArg1(:,:,x),AvgFil,"same");
    G_filtered_avg(:,:,x) = convn(inputArg2(:,:,x),AvgFil,"same");
    B_filtered_avg(:,:,x) = convn(inputArg3(:,:,x),AvgFil,"same");

end

%background light filtering

end

