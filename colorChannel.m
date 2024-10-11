ROI_R = zeros(1920,1080,30);
ROI_G = zeros(1920,1080,30);
ROI_B = zeros(1920,1080,30);

function [y] = myfun(data)
    
    for x = 1:30
        image1R = str(x);
        image1G = image1(:,:,2);
        image1B = image1(:,:,3);
        RegRowSt = 500;
        RegRowFin = 900;
        RegColSt = 500;
        RegColFin = 700;
        
        %getting the rgb colour channel
        image1GF = double(image1G);
        image1RF = double(image1R);
        image1BF = double(image1B);

end