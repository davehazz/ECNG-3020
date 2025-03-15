function [vectorRchannel,vectorGchannel,vectorBchannel,test_frames] = color_channel(vid)

%read vids
%this function reads in all the images into a vector
vid = VideoReader("IMG_9559.MOV");
len = 600;

LR=1920;
SR=640;
LC=1080;
SC= 480;

vectorRchannel = zeros(1080,1920,len);
vectorGchannel = zeros(1080,1920,len);
vectorBchannel = zeros(1080,1920,len);


test_frames = read(vid,[300 900]);


stopper = round(vid.Duration/5);
start = 1;
endd = round((1/stopper)*vid.NumFrames);

p=1;
%{
while stopper >=1
    test_frames = read(vid,[start endd]);
    for x =start:endd
    %this gave an error bacause i'm saving one image in a vector for size
    %30
    %str(x)= sprintf(" 'C:\\Users\\unity\\source\\repos\\Year two\\clas\\ECNG 3020 Heart Rate\\frame%d.jpg' ",count);
    %shell(x) = double(imread(str(x)));
    
    %vectorRchannel(:,:,x) = double(nz(:,:,x)); %aquire red channel

    vectorRchannel(:,:,start) = double(test_frames(:,:,1,start)); %aquire red channel
    vectorGchannel(:,:,start) = double(test_frames(:,:,2,start)); %aquire green channel
    vectorBchannel(:,:,start) = double(test_frames(:,:,3,start)); %aquire blue channel

    %channelG = nz(:,:,2);
    %shell(x) = double(nz);
    %shell(x) = double(read(v,x));
    
    end 
    fprintf('---Pass----%d',p);
    disp(mean(vectorRchannel));
    disp(mean(vectorGchannel));
    disp(mean(vectorBchannel));
    p=p+1;
    start = endd;
    endd = endd +endd;
    stopper = stopper-1;
end



%}

for x =1:600
    %this gave an error bacause i'm saving one image in a vector for size
    %30
    %str(x)= sprintf(" 'C:\\Users\\unity\\source\\repos\\Year two\\clas\\ECNG 3020 Heart Rate\\frame%d.jpg' ",count);
    %shell(x) = double(imread(str(x)));
    
    %vectorRchannel(:,:,x) = double(nz(:,:,x)); %aquire red channel

    vectorRchannel(:,:,x) = double(test_frames(:,:,1,x)); %aquire red channel
    vectorGchannel(:,:,x) = double(test_frames(:,:,2,x)); %aquire green channel
    vectorBchannel(:,:,x) = double(test_frames(:,:,3,x)); %aquire blue channel
    
  
    %channelG = nz(:,:,2);
    %shell(x) = double(nz);
    %shell(x) = double(read(v,x));

end 
%
%displaying the individual channels
%figure(1);imshow(uint8(imr));
%figure(2);imshow(uint8(img));
%figure(3);imshow(uint8(imb));
%showing the final picture
%figure(4);imshow(uint8(vectorRchannel(:,:,:,1)));

%OUTPUT FOR FUNCTION

%[ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel);


%{
im = vectorRchannel(:,:,1,1);
im = vectorRchannel(:,:,2,1);
img = vectorRchannel(:,:,2,1);
im = vectorRchannel(:,:,1,1);
imr = vectorRchannel(:,:,1,1);
img = vectorRchannel(:,:,2,1);
imb = vectorRchannel(:,:,3,1);
%}


end