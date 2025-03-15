%creating a video reader object
v = VideoReader("IMG_9559.MOV");

videor = readFrame(v);
imved = v.NumFrames;

te(:,:,1) = double(videor(:,:,1,1)); %aquire red channel
t = 0:1/29.98:20;
%-------------Part 1:vectors for each color channel---------

[vectorRchannel,vectorGchannel,vectorBchannel,imgg] = color_channel(v);

figure("Name","whole image");imshow(videor(:,:,2));

%figure(1);imshow(uint8(vectorRchannel(:,:,1)));
figure(1);imshow(uint8(vectorRchannel(:,:,1)));
figure(2);imshow(uint8(vectorGchannel(:,:,1)));
figure(3);imshow(uint8(vectorBchannel(:,:,1)));



%------------Part 2: Region of Interest--------------

[ROI_R,ROI_G,ROI_B] = region_of_interest(vectorRchannel,vectorGchannel,vectorBchannel);



[R_Source_m2,G_Source_m2,B_Source_m2] = signal_source(ROI_R,ROI_G,ROI_B);
R_m2_1 = zeros(1,600);
G_m2_2 = zeros(1,600);
B_m2_3 = zeros(1,600);

for x = 1:600
    R_m2_1(:,x) = (G_Source(:,x) ./ R_Source(:,x)) -1;
    G_m2_2(:,x) = (G_Source(:,x) ./ B_Source(:,x)) -1;
    B_m2_3(:,x) = (B_Source(:,x) ./ R_Source(:,x)) -1;

end

plot(t,R_m2_1);


figure(4);imshow(uint8(ROI_G(:,:,7)));
figure(5);imshow(uint8(ROI_G(:,:,60)));
figure(6);imshow(uint8(ROI_G(:,:,90)));


%------------Part 3:filtering the image----------------
[R_filtered,G_filtered,B_filtered] = filter_img_avg(ROI_R,ROI_G,ROI_B);

figure(7);imshow(uint8(G_filtered(:,:,10))); %at 10
figure(8);imshow(uint8(G_filtered(:,:,10)));% at 5
figure(9);imshow(uint8(G_filtered(:,:,10)));% at 0


%background light reduction
[ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B] = ROI_light_reduction(vectorRchannel,vectorGchannel,vectorBchannel);
figure(10);imshow(uint8(ROI_Light_filter_G(:,:,60)));

[RC_filtered,GC_filtered,BC_filtered] = filter_img_avg(ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B);


%------------Part 4: Signal Source------------------
[R_Source,G_Source,B_Source] = signal_source(R_filtered,G_filtered,B_filtered);
%[R_Source,G_Source,B_Source] = signal_source(ROI_R,ROI_G,ROI_B);

figure('Name','after bk reduction')
imshow(uint8(G_filtered(:,:,10)));

[R_Source_light,G_Source_light,B_Source_light] = signal_source(RC_filtered,GC_filtered,BC_filtered);
%without any moving average filter
%[R_Source_light,G_Source_light,B_Source_light] = signal_source(ROI_Light_filter_R,ROI_Light_filter_G,ROI_Light_filter_B);
j = medfilt2(RC_filtered,[3 3]);
i = median(RC_filtered,1);
%light reduction
R_c = zeros(1,600);
G_c = zeros(1,600);
B_c = zeros(1,600);


for x = 1:600
    R_c(:,x) = R_Source(:,x) - 0.75*R_Source_light(:,x);
    G_c(:,x) = G_Source(:,x) - 0.75*G_Source_light(:,x);
    B_c(:,x) = B_Source(:,x) - 0.75*B_Source_light(:,x);

end

figure(11);
plot(t,R_c,Color="r",Marker=".");
hold on
%figure(5);
plot(t,G_c,Color="g",Marker=".");
%figure(6);
plot(t,B_c,Color="b",Marker=".");
hold off


figure('Name','source signals');
plot(t,R_Source_light,Color="r",Marker=".");
hold on
%figure(5);
plot(t,G_Source_light,Color="g",Marker="^");
%figure(6);
plot(t,B_Source_light,Color="b",Marker="o");
hold off



%------------Part 5: Post Signal Filtering------------------
 
%Normalization
[normalizeR,normalizeG,normalizeB] = normalize_sig(R_c,G_c,B_c);

newR = normalize(R_Source);
newG = normalize(G_Source);
newB = normalize(B_Source);
figure(13);plot(t,newR);
%hold on
figure(14);plot(t,R_Source);
hold off


%before ICA
figure('Name','After Normalization');
plot(t,normalizeR,LineStyle="-");
hold on
plot(t,normalizeG,LineStyle="--");
plot(t,normalizeB,LineStyle=":");
hold off

%Detrending

R_de = zeros(1,700);
G_de = zeros(1,700);
B_de = zeros(1,700);

R_de = detrend(normalizeR);
G_de = detrend(normalizeG);
B_de = detrend(normalizeB);

%matlab normalize
RR_de = detrend(newR);
GG_de = detrend(newG);
BB_de = detrend(newB);

figure(16);
plot(t,R_de,LineStyle="-");
hold on
plot(t,RR_de,LineStyle=":");
legend('MINE  Data','MATLAB  Data')
hold off

figure(17);
plot(t,R_de);
hold on
plot(t,G_de);
plot(t,B_de);
hold off

%-------------more moving average filtering--------


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





%------------Part 6: ICA------------------

timmm = v.FrameRate*30;

%figure(2);plot(G_Source,LineStyle="--",Color="g");
%figure(3);plot(B_Source,LineStyle=":",Color="b");

combined_signal = ones(3,200);
combined_signal(1,:) = normalizeR(1,:);
combined_signal(2,:) = normalizeG(1,:);
combined_signal(3,:) = normalizeB(1,:);


tim = 0:0.03333:1;
%ica


Zica_G = fastICA(G_de,3,"kurtosis",0);%mine
Zica_GG = fastICA(GG_de,3,"kurtosis",0);%matlab
Zica_R = fastICA(R_de,3,"kurtosis",0);
Zica_RR = fastICA(RR_de,3,"kurtosis",0);
Zica_B = fastICA(B_de,3,"kurtosis",0);
ZICA_MAF_R = fastICA(y,3,"kurtosis,0",0);
ZICA_MAF_G = fastICA(y1,3,"kurtosis,0",0);


%Zica_RR = fastICA(RR_de,3);

figure(19);
plot(t,Zica_R(1,:),Color="r");
hold on
plot(t,Zica_G(1,:),Color='g');
plot(t,Zica_B(1,:),Color='b');
hold off


figure(20);
plot(t,Zica_R(1,:),LineStyle="-",Color="r");
hold on
plot(t,Zica_G(1,:),LineStyle=":",Color='g');
plot(t,Zica_B(3,:),LineStyle="-.",Color='b');


figure(21);
plot(t,Zica_G(1,:),LineStyle="-",Color='r');
hold on
plot(t,Zica_G(2,:),LineStyle=":",Color='g');
plot(t,Zica_G(3,:),LineStyle="-.",Color='b');


figure(22);
plot(t,Zica_B(1,:),LineStyle="-",Color='r');
hold on
plot(t,Zica_B(2,:),LineStyle=":",Color='g');
plot(t,Zica_B(3,:),LineStyle="-.",Color='b');



%------------Part 7: FFT and PSD------------------

figure(3);pwelch(Zica_R(1,:));
figure(1);pwelch(Zica_G(1,:));
figure(2);pwelch(R_Source);
figure(4);pwelch(normalizeG);
[ppx,f] = pwelch(Zica_G(1,:));
pwelch(normalizeG);


%----------FFT OF the signal-------------
Fs = 1000;            % Sampling frequency                    
T = 0.001;             % Sampling period       
L = 100;             % Length of signal
time_vec = (0:L-1);        % Time vector
new_time_vec =  time_vec.*T;


Y = fft(Zica_G(1,:));
plot(Fs/L*(0:L-1),abs(Y),"LineWidth",2)
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")



P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs/L*(0:(L/2));
plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")


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


ft1 = A*sin(w*ts);
ft2 = A*sin(2*pi*300*ts);
ft3 = A*sin(2*pi*150*ts);
ft4 = A*sin(2*pi*100*ts);

plot(ts,ft1);
hold on
plot(ts,ft2);
plot(ts,ft3);
plot(ts,ft4);

hold off

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

Fs = 1000; % Sampling frequency (Hz)
Fc1 = 0.8; % Lower cutoff frequency (Hz)
Fc2 = 2; % Upper cutoff frequency (Hz)
[b, a] = butter(4, [Fc1 Fc2]/(Fs/2)); % 4th order Butterworth filter

% Apply the filter to a signal
filtered_signal = filter(b, a, Zica_G(1,:));
plot(filtered_signal);

figure(1);bandpass(ZICA_MAF_G(1,:), [0.8 2], 300);%
figure(2);bandpass(ZICA_MAF_R(1,:), [0.8 2], 300);%
figure(30);bandpass(Zica_R(1,:), [0.8 2], 100);
figure(40);bandpass(Zica_G(1,:), [0.8 2], 300);
figure(6);bandpass(Zica_G(2,:), [0.8 2], 500);



hold off

%----------PSD Of the signal---------



f = Fs*(0:200/2-1)/200;
S_meg = abs(s_oneSide)/(100);

plot(S_meg);
plot(Zica_B(1,:));


ppx = pwelch(S_meg);
plot(ppx); 

y = bandpass(ppx,[1 2], 1000);
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
%t = 0:2*pi/1000:2*pi;

Fs= 1000; %sampling frequency
T = 1/Fs; %sampling period
L = 3000;
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
    seg(x,:)= y(m:n);
    m=n;
    n=n+n;
end

%plotting the full time period
figure("Name","full time");
plot(t,y);

%finding the pwd using the fft
nfft = fft(seg_1);%the fft

plot(Fs/L*(0:500-1),abs(nfft),"LineWidth",3)
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")

real_val = abs(nfft);%magnitudes of it

%length of it
length_nfft = length(abs(nfft));
sq_val = (abs(nfft).^2)./(2);  %the actual PSD
%sq_val = (abs(nfft).^2)./(length(nfft));  %the actual PSD

figure("Name","seg 1 psd");
plot(sq_val);

%{
    % i am looking for the median peak within the 0.8 to 2Hz range

%}
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

IFaces = insertObjectAnnotation(I,'rectangle',bboxes,'Face');   
figure
imshow(IFaces)
title('Detected faces');