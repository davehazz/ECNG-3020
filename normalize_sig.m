function [outputArg1,outputArg2,outputArg3] = normalize_sig(inputArg1,inputArg2,inputArg3)
%UNTITLED Summary of this function goes here
%   this function provides the normilization of the source signal
% it normalize = ((Sx - mean) / standard deviation )

c = 600;

outputArg1 = zeros(1,c);
outputArg2 = zeros(1,c);
outputArg3 = zeros(1,c);

for start = 1:c
    outputArg1(start) = ( (inputArg1(start) - mean(inputArg1)) / std(inputArg1));
    outputArg2(start) = ( (inputArg2(start) - mean(inputArg2)) / std(inputArg2));
    outputArg3(start) = ( (inputArg3(start) - mean(inputArg3)) / std(inputArg3));

end