%creating a video reader object
v = VideoReader("IMG_9136.MOV");
videor = readFrame(v);
imved = v.NumFrames;

te(:,:,1) = double(videor(:,:,1,1)); %aquire red channel
t = 0:1/29.97:5.3717;
%-------------Part 1:vectors for each color channel---------

[vectorRchannel,vectorGchannel,vectorBchannel,imgg] = color_channel(v);

%figure(1);imshow(uint8(vectorRchannel(:,:,1)));
figure(1);imshow(uint8(vectorGchannel(:,:,1)));
%figure(3);imshow(uint8(vectorBchannel(:,:,1)));


%------------Part 2: Region of Interest--------------

[ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel);

%[ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B] = ROI_light_reduction(vectorRchannel,vectorGchannel,vectorBchannel);

figure(2);imshow(uint8(ROI_G(:,:,1)));


%------------Part 3:filtering the image----------------
[R_filtered,G_filtered,B_filtered] = filter_img_avg(ROI_R,ROI_G,ROI_B);

figure(3);imshow(uint8(G_filtered(:,:,1)));


%------------Part 4: Signal Source------------------
[R_Source,G_Source,B_Source] = signal_source(R_filtered,G_filtered,B_filtered);


for x = 1:161
    R_Source(x) = mean2(inputArg1(:,:,x));
    G_Source(x) = mean2(inputArg1(:,:,x));
    B_Source(x) = mean2(inputArg1(:,:,x));

end

figure(4);plot(t,R_Source,Color="r");
figure(5);plot(t,G_Source,Color="g");
figure(6);plot(t,B_Source,Color="b");



%------------Part 5: Post Signal Filtering------------------
 
%Normalization
[normalizeR,normalizeG,normalizeB] = normalize_sig(R_c,G_c,B_c);

newR = normalize(R_Source);
figure(1);plot(t,newR);
%hold on
figure(2);plot(t,R_Source);
hold off


%before ICA
plot(t,normalizeR,LineStyle="-");
hold on
figure(7);plot(t,normalizeG,LineStyle="--");
plot(t,normalizeB,LineStyle=":");
hold off

%Detrending

R_de = zeros(1,700);
G_de = zeros(1,700);
B_de = zeros(1,700);

R_de = detrend(normalizeR);
G_de = detrend(normalizeG);
B_de(:,x) = detrend(normalizeB);

plot(t,R_de,LineStyle="-");
hold on
plot(t,normalizeR,LineStyle="-");
legend('Original Data','Detrended Data')

hold off


%------------Part 6: ICA------------------

timmm = v.FrameRate*30;

%figure(2);plot(G_Source,LineStyle="--",Color="g");
%figure(3);plot(B_Source,LineStyle=":",Color="b");

combined_signal = ones(3,200);
combined_signal(1,:) = normalizeR(1,:);
combined_signal(2,:) = normalizeG(1,:);
combined_signal(3,:) = normalizeB(1,:);


tim = 0:0.03333:6.66;
%ica


Zica_R = fastICA(normalizeR,3,"kurtosis",1);
Zica_G = fastICA(normalizeG,3,"kurtosis",0);
Zica_B = fastICA(normalizeB,3,"kurtosis",0);


figure(1);plot(t,Zica_R(1,:),LineStyle="-",Color="r");
hold on
plot(t,Zica_G(1,:),LineStyle="--",Color='g');
plot(t,Zica_B(1,:),LineStyle=":",Color='b');

figure(2);plot(tim,Zica_R(2,:),LineStyle="-",Marker="o",Color="r");
hold on
plot(tim,Zica_G(2,:),LineStyle=":",Marker="v",Color='g');
plot(tim,Zica_B(2,:),LineStyle="-.",Marker="*",Color='b');


figure(3);plot(tim,Zica_R(3,:),LineStyle="-",Marker="o",Color="r");
hold on
plot(tim,Zica_G(3,:),LineStyle=":",Marker="v",Color='g');
plot(tim,Zica_B(3,:),LineStyle="-.",Marker="*",Color='b');




%------------Part 7: FFT and PSD------------------

figure(3);pwelch(Zica_G(1,:));
figure(1);pwelch(G_Source);
figure(2);pwelch(R_Source);
figure(4);pwelch(normalizeG);
pwelch(normalizeG);

ff_R = fft2(Zica_G(1,:));
ff_g = fft2(normalizeG);
plot(abs(ff_R));

%----------FFT OF the signal-------------
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 500;             % Length of signal
T = (0:L-1)*T;        % Time vector

Y = fft(Zica_G(1,:));
plot(Fs/L*(0:L-1),abs(Y),"LineWidth",2)
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")

%test
A = 5;
f = 1/pi; %this is the frequency so the period will be 1/f
theta = 0; %this is the phase shift
ts = linspace(0,100,10000);  
%linspace goes from 0 to pi with 100000 subdivisions
w = 2*pi*f; %this is the angluar frequency or the fundamental frequency


ft1 = A*sin(w*ts);
ft2 = A*sin(2*pi*300*ts);
ft3 = A*sin(2*pi*150*ts);
ft4 = A*sin(2*pi*100*ts);

plot(ts,ft1);
hold on
plot(ts,ft2);
plot(ts,ft3);
plot(ts,ft4);

nfft = length(Zica_G(1,:));
periodogram(Zica_G(1,:),nfft,);

X = ft1+ft2+ft3+ft4;
plot(ts,X);

%showing fft
plot(ts,ft1);
Fs= 1000; %sampling frequency
T = 1/Fs; %sampling period
L = 10000;
t = (0:L-1)*T; %time vector

ft_new = fft(X);
plot(Fs/L*(0:L-1),abs(ft_new));

%----------Bandpass filtering--------

figure(1);bandpass(X, [50 200], 1000);
figure(2);bandpass(Zica_G(1,:), [1 2], 1000);
figure(3);bandpass(Zica_B(1,:), [1 2], 1000);
figure(4);bandpass(normalizeG, [1 2], 1000);

hold off

%----------PSD Of the signal---------

bandpass(abs(s_oneSide),[1 4],1000);

Fs = 10000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 6.67;             % Length of signal
t = (0:L-1)*T;        % Time vector

f = Fs*(0:200/2-1)/200;
S_meg = abs(s_oneSide)/(100);

plot(S_meg);
plot(Zica_B(1,:));


ppx = pwelch(S_meg);
plot(ppx); 

y = bandpass(ppx,[0.8 2], 1000);
plot(y);


y = fft(R_Source);
plot(t,abs(y));
figure(4);plot(y);

%{
%{
for x = 1:10
    %fprintf('hello\n');
    %sets each frame into the respective array
    str(x)= sprintf(" 'C:\\Users\\unity\\source\\repos\\Year two\\clas\\ECNG 3020 Heart Rate\\frame%d.jpg' ",count);
    %disp(str(x));
    %image_vec(x) = imread(str(x)); 
    count= count+1;

end 

%testim = imread(str(1));
%}






%------------PART 2: RGB Channel------------------ 

location = 'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame0.jpg';
image1 = imread(location);
image1R = image1(:,:,1);
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


%-------------Part 3: Select a ROI---------------------

%selecting the region of interest
ROI = image1GF(RegRowSt:RegRowFin,RegColSt:RegColFin);
ROI2 = image1RF(RegRowSt:RegRowFin,RegColSt:RegColFin);
ROI3 = image1BF(RegRowSt:RegRowFin,RegColSt:RegColFin);


%displaying the region of interest
figure(1);imshow(uint8(ROI));
figure(2);imshow(uint8(ROI2));
figure(3);imshow(uint8(ROI3));
%}

%-------------Test Part---------------
%{
location = 'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame0.jpg';
image1G = image1d(:,:,2);
RegRowSt = 500;
RegRowFin = 900;
RegColSt = 500;
RegColFin = 700;
image1GF = double(image1G);
ROI = image1GF(RegRowSt:RegRowFin,RegColSt:RegColFin);
%}


%---------------PART 4: Use Filter---------------
%{
%image_test =  'C:\Users\unity\source\repos\Year two\clas\ECNG 3020 Heart Rate\frame9.jpg';
%image1 = imread(RO1);
%image1d = double(image1);


%%% Averaging Filter
%Task 3: Create a 2D averaging filter of size n
%create a variable with the filter size
AvgFilSize = 50;
%Create a 2D empty shell with ones with the size
AvgFil = ones(AvgFilSize,AvgFilSize);
% Make sure the filter sums up to one (normalization)
AvgFil = AvgFil/sum(AvgFil(:));
% apply the filter to the image using convolution
ImgR = convn(double(ROI2),AvgFil,'same');
ImgG = convn(double(ROI),AvgFil,'same');
ImgB = convn(double(ROI3),AvgFil,'same');

% show the image
figure(7);subplot(2,2,3);imshow(uint8(ImgF)); title('ImgF-10');
figure(7);subplot(2,2,4);imshow(uint8(ROI)); title('-original');

%}

%-----------PART 5: Independent Component Analysis---------------

%-----------PART 5, a)signal source creation---------------
% 
%i am finding the mean of the image
% i am doing this to obtain a time series of means of the frames
%{
meanIntensityValueR = mean2(ImgR);
meanIntensityValueG = mean2(ImgG);
meanIntensityValueB = mean2(ImgB);
 %}


%{
%comp = 2;
%train ICA model
y = fft2(image1GF);
y = abs(y);
figure(2);
y = log(y+1);

imshow(y,[]);
title('fourier transform of an image')

Mld = rica(y,2);

%plot(Mld.)
%apply ICA
%data_ICA = transform(Mld,)
%}



%}
