function [R_filtered_light,G_filtered_light,B_filtered_light] = filter_img_light()
%UNTITLED4 Summary of this function goes here
%this is for removing background light

%   this removes background light


%this function performs the filering on the image
R_filtered_light = zeros(401,201,200);
G_filtered_light = zeros(401,201,200);
B_filtered_light = zeros(401,201,200);

[R_Source,G_Source,B_Source] = signal_source(R_filtered_avg,G_filtered_avg,B_filtered_avg);


for x = 1:200
    
    R_filtered_light(x) = convn(inputArg1(:,:,x),AvgFil,'same');
    G_filtered_light(x) = convn(inputArg1(:,:,x),AvgFil,'same');
    B_filtered_light(x) = convn(inputArg1(:,:,x),AvgFil,'same');

end


end