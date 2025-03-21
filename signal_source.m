function [R_Source,G_Source,B_Source] = signal_source(inputArg1,inputArg2,inputArg3)
%UNTITLED3 Summary of this function goes here

%   Detailed explanation goes here
%this function creates the signal source from finding the mean pixel
%intensities of each frame

size = 600;

R_Source = zeros(1,size);
G_Source = zeros(1,size);
B_Source = zeros(1,size);

for x = 1:size
    R_Source(x) = mean(inputArg1(:,:,x),"all");
    G_Source(x) = mean(inputArg1(:,:,x),"all");
    B_Source(x) = mean(inputArg1(:,:,x),"all");

end