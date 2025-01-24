function [ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B] = ROI_light_reduction(inputArg1,inputArg2,inputArg3)
%UNTITLED Summary of this function goes here
%   this selects the region of interest for the background light reduction


a = 101;
b = 101;
c = 161;

ROI_Light_filter_R = zeros(a,b,c);
ROI_Light_filter_G = zeros(a,b,c);
ROI_Light_filter_B = zeros(a,b,c);

        RegRowSt = 800;
        RegRowFin = 900;
        RegColSt = 500;
        RegColFin = 600;


    for x = 1:161
        ROI_Light_filter_R(:,:,x) = inputArg1(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_Light_filter_G(:,:,x) = inputArg2(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_Light_filter_B(:,:,x) = inputArg3(RegRowSt:RegRowFin,RegColSt:RegColFin,x);

        %getting the rgb colour channel

    end




end