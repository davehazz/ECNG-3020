function [ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel)
a = 41;
b = 31;
c = 600;

ROI_R = zeros(a,b,c);
ROI_G = zeros(a,b,c);
ROI_B = zeros(a,b,c);

        RegRowSt = 550;
        RegRowFin = 590;
        RegColSt = 1030;
        RegColFin = 1060;


    for x = 1:600
        ROI_R(:,:,x) = vectorRchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_G(:,:,x) = vectorGchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_B(:,:,x) = vectorBchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);

        %getting the rgb colour channel

    end
%figure(2);imshow(uint8(ROI_G(:,:,1)));

end