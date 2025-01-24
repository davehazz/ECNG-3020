function [R_Source,G_Source,B_Source] = signal_source(inputArg1,inputArg2,inputArg3)
%UNTITLED3 Summary of this function goes here

%   Detailed explanation goes here
%this function creates the signal source from finding the mean pixel
%intensities of each frame

R_Source = zeros(1,161);
G_Source = zeros(1,161);
B_Source = zeros(1,161);

for x = 1:161
    R_Source(x) = mean2(inputArg1(:,:,x));
    G_Source(x) = mean2(inputArg1(:,:,x));
    B_Source(x) = mean2(inputArg1(:,:,x));

end