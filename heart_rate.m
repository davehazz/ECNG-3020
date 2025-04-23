%creating a video reader object
v = VideoReader("C:\Users\unity\Documents\Engineering uWI\YEAR THREE\ECNG 3020\project code\videos\IMG_9988.MOV");


videor = read(v, [1 600]);
imved = v.NumFrames;


te(:,:,1) = double(videor(:,:,1,1)); %aquire red channel
t = 0:1/29.98:20;
%-------------Part 1:vectors for each color channel---------

%-- forget this --> [vectorRchannel,vectorGchannel,vectorBchannel,imgg] = color_channel(v);

%--- run from here 
len = 600;
vectorRchannel = zeros(1088,1920,len);
vectorGchannel = zeros(1088,1920,len);
vectorBchannel = zeros(1088,1920,len);
%--- run to here 

currAxes = axes;
vidFrame = readFrame(v);
image(vidFrame,"Parent",currAxes);


%--- run this line
for x =1:600
    
    vectorRchannel(:,:,x) = double(videor(:,:,1,x)); %aquire red channel
    vectorGchannel(:,:,x) = double(videor(:,:,2,x)); %aquire green channel
    vectorBchannel(:,:,x) = double(videor(:,:,3,x)); %aquire blue channel

end 
%--- run this line


figure("Name","whole image");
imshow(videor(:,:,2));

%figure(1);imshow(uint8(vectorRchannel(:,:,1)));
figure('Name','red channel');imshow(uint8(vectorRchannel(:,:,1)));
figure('Name','green channel');imshow(uint8(vectorGchannel(:,:,1)));
figure('Name','blue chennel');imshow(uint8(vectorBchannel(:,:,1)));



%------------Part 2: Region of Interest--------------

vectorRchannel = 0;
vectorGchannel = 0;
vectorBchannel = 0;
[ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel); %--- run this line

figure(4);imshow(uint8(ROI_G(:,:,10))); %--- run this line


[R_Source_m2,G_Source_m2,B_Source_m2] = signal_source(ROI_R,ROI_G,ROI_B); %--- run this line

%{
R_m2_1 = zeros(1,600);
G_m2_2 = zeros(1,600);
B_m2_3 = zeros(1,600);

for x = 1:600
    R_m2_1(:,x) = (G_Source_m2(:,x) ./ R_Source_m2(:,x)) -1;
    G_m2_2(:,x) = (G_Source_m2(:,x) ./ B_Source_m2(:,x)) -1;
    B_m2_3(:,x) = (B_Source_m2(:,x) ./ R_Source_m2(:,x)) -1;

end

figure('Name','ratio method');
plot(t,R_m2_1,Color='r');
hold on
plot(t,G_m2_2,Color='g');
plot(t,B_m2_3,Color='b');
hold off
%}


R_m2_1 = 0;
G_m2_2 = 0;
B_m2_3 = 0;


figure(5);imshow(uint8(ROI_G(:,:,60)));
figure(6);imshow(uint8(ROI_G(:,:,90)));


brk = uint8(ROI_G(:,:,10));
brk = brk./3;
imshow(brk);

ROI_R = 0;
ROI_G = 0;
ROI_B = 0;

%------------Part 3:filtering the image----------------
[R_filtered,G_filtered,B_filtered] = filter_img_avg(ROI_R,ROI_G,ROI_B); %--- run this line

figure(7);imshow(uint8(G_filtered(:,:,10))); %at 10
figure(8);imshow(uint8(G_filtered(:,:,10)));% at 5
figure(9);imshow(uint8(G_filtered(:,:,10)));% at 0


%background light reduction
[ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B] = ROI_light_reduction(vectorRchannel,vectorGchannel,vectorBchannel); %--- run this line   
figure(10);imshow(uint8(ROI_Light_filter_G(:,:,10)));

[RC_filtered,GC_filtered,BC_filtered] = filter_img_avg(ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B);%--- run this line
figure(9);imshow(uint8(GC_filtered(:,:,10)));% at 0


%------------Part 4: Signal Source------------------
[R_Source,G_Source,B_Source] = signal_source(R_filtered,G_filtered,B_filtered);%--- run this line
%[R_Source,G_Source,B_Source] = signal_source(ROI_R,ROI_G,ROI_B);

figure('Name','after bk reduction')
imshow(uint8(G_filtered(:,:,90)));

[R_Source_light,G_Source_light,B_Source_light] = signal_source(RC_filtered,GC_filtered,BC_filtered);%--- run this line
%without any moving average filter
%[R_Source_light,G_Source_light,B_Source_light] = signal_source(ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B);
j = medfilt2(RC_filtered,[3 3]);
i = median(RC_filtered,1);

%----run below lines
R_c = zeros(1,600);
G_c = zeros(1,600);
B_c = zeros(1,600);
%--- run this line


%--- run this line
for x = 1:600
    R_c(:,x) = R_Source(:,x) - 0.75*R_Source_light(:,x);
    G_c(:,x) = G_Source(:,x) - 0.75*G_Source_light(:,x);
    B_c(:,x) = B_Source(:,x) - 0.75*B_Source_light(:,x);

end
%---------till here

%{
figure(11);
plot(t,R_c,Color="r");
hold on
%figure(5);
plot(t,G_c,Color="g");
%figure(6);
plot(t,B_c,Color="b");
hold off
%}

%---------MEAN vs MEDIAN Pixel Intensity--------

%{
G_c_mean = G_c;
G_c_median;

plot(t,G_c_mean);
hold on
plot(t,G_c_median);
legend('Mean  Data','Median  Data')
hold off



figure('Name','source signals');
plot(t,R_Source_light,Color="r",Marker=".");
hold on
%figure(5);
plot(t,G_Source_light,Color="g",Marker="^");
%figure(6);
plot(t,B_Source_light,Color="b",Marker="o");
hold off

%}




%------------Part 5: Post Signal Filtering------------------
 
%Normalization
[normalizeR,normalizeG,normalizeB] = normalize_sig(R_c,G_c,B_c); %--- run this line

%{
R_c = 0;
G_c = 0;
B_c = 0;

newR = normalize(R_Source);
newG = normalize(G_Source);
newB = normalize(B_Source);
%my normalize was not within the range of -2 to +2
figure('Name','after normalization');
plot(t,newG);
hold on
plot(t,normalizeG);
legend('Matlab Data','My Data')
hold off



%before ICA
figure('Name','After Normalization');
plot(t,normalizeR,Color="r");
hold on
figure('Name','After Normalization');
plot(t,normalizeG,Color='g');
plot(t,normalizeB,Color="b");
hold off

figure('Name','before detreding');
plot(t,normalizeG,LineStyle="-");
hold on
plot(t,newG,LineStyle=":");
legend('MINE  Data','MATLAB  Data')
hold off

%}




%Detrending


%-------run below lines
R_de = detrend(normalizeR);
G_de = detrend(normalizeG);
B_de = detrend(normalizeB);
%-------till here


%matlab normalize
%{

RR_de = detrend(newR);
GG_de = detrend(newG);
BB_de = detrend(newB);

RR_de = 0;
GG_de = 0;
BB_de = 0;


figure('Name','detrending');
plot(t,G_de,LineStyle="-");
hold on
plot(t,GG_de,LineStyle=":");
legend('MINE  Data','MATLAB  Data')
hold off

figure(17);
plot(t,R_de);
hold on
plot(t,G_de);
plot(t,B_de);
hold off

%}

%-------------more moving average filtering--------


%------------- STEP NOT NEEDED-------------
%{
windowSize = 7; %setting filter settings
b = (1/windowSize)*ones(1,windowSize);
a = 1;

%applying the filer
y = filter(b,a,RR_de);
y1 = filter(b,a,GG_de);


figure(18);
plot(t,y,LineStyle="-");
hold on
plot(t,y1,LineStyle="-");
legend('R Channel','G channel')
hold off


%}



%------------Part 6: ICA------------------

%ica

Zica_G = fastICA(G_de,3,"kurtosis",0); %--------run this line

%{
Zica_GG = fastICA(G_de,5,"negentropy",0);%matlab
Zica_R = fastICA(R_de,3,"kurtosis",0);
Zica_RR = fastICA(RR_de,3,"kurtosis",0);
Zica_B = fastICA(B_de,3,"kurtosis",0);
ZICA_BB = fastICA(BB_de,3,"kurtosis,0",0);



%Zica_RR = fastICA(RR_de,3);

figure(19);
plot(t,Zica_GG(1,:),Color="r");
hold on
plot(t,Zica_G(2,:),Color='g');
plot(t,Zica_GG(5,:),Color='b');
hold off


figure('Name','After ICA');
plot(t,Zica_R(1,:),LineStyle="-",Color="r");
hold on
plot(t,Zica_G(1,:),Color='g');
plot(t,Zica_B(1,:),LineStyle="-.",Color='b');


figure(21);
plot(t,Zica_G(1,:),Color='g',LineStyle='-');
hold on
plot(t,Zica_G(3,:),Color='r',LineStyle='--');
plot(t,Zica_G(5,:),Color='g',LineStyle=':');


figure(22);
plot(t,Zica_GG(3,:),LineStyle="-",Color='r');
hold on
plot(t,Zica_GG(3,:),LineStyle=":",Color='g');
plot(t,Zica_GG(3,:),LineStyle="-.",Color='b');


te = snr(Zica_G(1,:));
te = snr(Zica_G(2,:));
%}



%------------Part 7: FFT and PSD------------------

%{
nnn = bandpass(Zica_G(1,:),[0.8 2],200);
figure(3);pwelch(Zica_R(1,:));
figure(1);pwelch(Zica_G(1,:)); 
figure(2);pwelch(R_Source);


figure(4);pwelch(normalizeG);

%}

%----run this line from here
[ppxg,fG] = pwelch(nnn);
[ppxgg,fGG] = pwelch(Zica_G(1,:));

figure('Name','signal');
plot(fGG,10*log10(ppxgg));
%--------- to there

pwelch(normalizeG);


%----------FFT OF the signal-------------

%{
N =1000;
Y = Y(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(Y).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:2*pi/N:pi;

plot(freq/pi,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Normalized Frequency (\times\pi rad/sample)")
ylabel("Power/Frequency (dB/(rad/sample))")
%test
A = 5;
f = 1/pi; %this is the frequency so the period will be 1/f
theta = 0; %this is the phase shift
ts = linspace(0,100,10000);  
%linspace goes from 0 to pi with 100000 subdivisions
w = 2*pi*f; %this is the angluar frequency or the fundamental frequency

%}



%---------------------Bandpass filtering------------------


figure(40);bandpass(Zica_G(1,:), [0.8 2], 200); %----run this line

%{
figure(30);bandpass(Zica_G(1,:), [0.8 2], 100);
figure(46);bandpass(Zica_GG(1,:), [0.8 2], 200);
figure(6);bandpass(Zica_G(2,:), [0.8 2], 500);
%}


%{



%t = 0:2*pi/1000:2*pi;

Fs= 1000; %sampling frequency
T = 1/Fs; %sampling period
L = 20000;
t = (0:L-1)*T; %time vector

y = 0.8 + 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

%plotting the forst 0.5 seconds
figure("Name","first 0.5 seconds");
plot(t(1:500),y(1:500));

seg_1 = y(1:500);
seg = zeros(6,500);
n = 500;
m=1;

for x=1:6
    seg(x,:)= y(x,m:n);
    m=n;
    n=n+n;
end

%plotting the full time period
figure("Name","full time");
plot(t,y);

%finding the pwd using the fft
nnn = bandpass(Zica_G(1,:),[0.8 2],200);
nfft_1 = fft(nnn(1:150));%the fft
nfft_2 = fft(nnn(151:300));%the fft
nfft_3 = fft(nnn(301:450));%the fft
nfft_4 = fft(nnn(451:600));%the fft




plot(Fs/L*(0:600-1),abs(nfft))
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")


%length of it

sq_val_1 = (abs(nfft_1).^2)./(length(nfft_1));  %the actual PSD
sq_val_2 = (abs(nfft_2).^2)./(length(nfft_2));  %the actual PSD
sq_val_3 = (abs(nfft_3).^2)./(length(nfft_3));  %the actual PSD
sq_val_4 = (abs(nfft_4).^2)./(length(nfft_4));  %the actual PSD



figure("Name","seg 1 psd");
plot(Fs/L*(0:150-1),sq_val_1);
hold on
plot(Fs/L*(151:300),sq_val_2);
plot(Fs/L*(301:450),sq_val_3);
plot(Fs/L*(451:600),sq_val_4);
hold off

sq = zeros(4,150);


sq(1,:) = sq_val_1;
sq(2,:) = sq_val_2;
sq(3,:) = sq_val_3;
sq(4,:) = sq_val_4;

fin = median(sq);

plot(Fs/L*(0:150-1),fin);






Y = fft(Zica_G(1,:));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")


bfd = bandpass(seg_1,[0.8 2],100);
plot(t(1:500),bfd);
st =  fft(bfd);
mag = abs(st).^2;

plot(Fs/L*(0:250-1),mag);


ppx = periodogram(y(1:500),[],50,"onesided");

plot(Fs/L*(0:L-1)/.2,abs(nfft));

faceDetector = vision.CascadeObjectDetector;
I = imread("C:\Users\unity\Documents\Engineering uWI\YEAR THREE\ECNG 3020\project code\face.jpg");
bboxes = faceDetector(I);


% Read a video frame and run the face detector.
bbox = step(faceDetector, I);

% Draw the returned bounding box around the detected face.
I = insertShape(I, "rectangle", bbox);
figure; imshow(I); title("Detected face");

bboxPoints = bbox2points(bbox(1, :));

points = detectMinEigenFeatures(im2gray(I),"ROI",  bbox);

% Display the detected points.
figure, imshow(videoFrame), hold on, title("Detected features");
plot(points);

%}
