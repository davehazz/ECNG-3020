function [ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel)
a = 151;
b = 301;
c = 700;

ROI_R = zeros(a,b,c);
ROI_G = zeros(a,b,c);
ROI_B = zeros(a,b,c);

        RegRowSt = 500;
        RegRowFin = 650;
        RegColSt = 300;
        RegColFin = 600;


    for x = 1:700
        ROI_R(:,:,x) = vectorRchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_G(:,:,x) = vectorGchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_B(:,:,x) = vectorBchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);

        %getting the rgb colour channel

    end
%figure(2);imshow(uint8(ROI_G(:,:,1)));

end