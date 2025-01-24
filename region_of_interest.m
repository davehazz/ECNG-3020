function [ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel)
a = 151;
b = 401;
c = 161;

ROI_R = zeros(a,b,c);
ROI_G = zeros(a,b,c);
ROI_B = zeros(a,b,c);

        RegRowSt = 450;
        RegRowFin = 600;
        RegColSt = 400;
        RegColFin = 800;


    for x = 1:161
        ROI_R(:,:,x) = vectorRchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_G(:,:,x) = vectorGchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_B(:,:,x) = vectorBchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);

        %getting the rgb colour channel

    end

end