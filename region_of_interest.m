function [ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel)
  
ROI_R = zeros(401,201,200);
ROI_G = zeros(401,201,200);
ROI_B = zeros(401,201,200);

        RegRowSt = 500;
        RegRowFin = 900;
        RegColSt = 500;
        RegColFin = 700;


    for x = 1:200
        ROI_R(:,:,x) = vectorRchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_G(:,:,x) = vectorGchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);
        ROI_B(:,:,x) = vectorBchannel(RegRowSt:RegRowFin,RegColSt:RegColFin,x);

        %getting the rgb colour channel

    end

end