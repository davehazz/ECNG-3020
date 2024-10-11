%read vids
%this function reads in all the images into a vector

vectorRchannel = zeros(1920,1080,3,30);
vectorGchannel = zeros(1920,1080,30,"double");
vectorBchannel = zeros(1920,1080,30,"double");

%iuiy = cast(video,double);


nelocation = 'C:\\Users\\unity\\source\\repos\\Year two\\clas\\ECNG 3020 Heart Rate\\frame0.jpg';

vv = VideoReader("IMG_6010.MOV");
test_frames = read(vv,[1 30]);

for x =1:30
    %this gave an error bacause i'm saving one image in a vector for size
    %30
    %str(x)= sprintf(" 'C:\\Users\\unity\\source\\repos\\Year two\\clas\\ECNG 3020 Heart Rate\\frame%d.jpg' ",count);
    %shell(x) = double(imread(str(x)));
    

    %vectorRchannel(x) = double(nz(:,:,1)); %aquire red channel

    vectorRchannel(:,:,:,x) = double(test_frames(:,:,:,x)); %aquire red channel
    %vectorGchannel(:,:,x) = double(nz(:,:,2)); %aquire green channel
    %vectorBchannel(:,:,x) = double(nz(:,:,3)); %aquire blue channel

    %channelG = nz(:,:,2);
    %shell(x) = double(nz);
    %shell(x) = double(read(v,x));

end 

%displaying the individual channels
figure(1);imshow(uint8(imr));
figure(2);imshow(uint8(img));
figure(3);imshow(uint8(imb));
%showing the final picture
figure(4);imshow(uint8(vectorRchannel(:,:,:,1)));


