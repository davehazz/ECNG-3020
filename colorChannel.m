ROI_R = zeros(1920,1080,30);
ROI_G = zeros(1920,1080,30);
ROI_B = zeros(1920,1080,30);


function [ROI_R,ROI_G,ROI_B] = myfun(vectorRchannel,vectorGchannel,vectorBchannel)
   
        RegRowSt = 500;
        RegRowFin = 900;
        RegColSt = 500;
        RegColFin = 700;


    for x = 1:30
        ROI_R(x) = vectorRchannel(x(RegRowSt:RegRowFin,RegColSt:RegColFin));
        ROI_G(x) = vectorRchannel(x(RegRowSt:RegRowFin,RegColSt:RegColFin));
        ROI_B(x) = vectorRchannel(x(RegRowSt:RegRowFin,RegColSt:RegColFin));

        %getting the rgb colour channel

    end

end