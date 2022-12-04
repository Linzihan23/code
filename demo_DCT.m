% =========================================================================
%Implementation of image noise reduction using sparse dictionaries
%Environment: Win7, Matlab2020a
%Time: 2022-12-3
% =========================================================================

%% Data reading and pre-processing
load('date.mat')
bb=8;
RR=4;
K=RR*bb^2; 
sigma = 25;

newrealphase=zeros(512);
for i=1:1:512
    x=realphase(:,i);
    t=1:1:512;    
    snr=1;
    px_dBW=1;
    y1=awgn(x,snr,px_dBW);
    newrealphase(:,i)=y1;
end

PSNRIn = 20*log10(255/sqrt(mean((newrealphase(:)-realphase(:)).^2)));

%% Gaussian random noise
% % newrealphase=realphase+sigma*randn(size(realphase));
% % PSNRIn = 20*log10(255/sqrt(mean((newrealphase(:)-realphase(:)).^2)));

%% DCT Denoising Test
[IoutDCT,output] = denoiseImageDCT(newrealphase, sigma, K);
PSNROut = 20*log10(255/sqrt(mean((IoutDCT(:)-realphase(:)).^2)));

figure;
subplot(1,2,1); imshow(newrealphase,[]); title(strcat(['Noisy image, ',num2str(PSNRIn),'dB']));
subplot(1,2,2); imshow(IoutDCT,[]); title(strcat(['Clean Image by DCT dictionary, ',num2str(PSNROut),'dB']));

figure;
I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb,0);
title('The DCT dictionary');

%% Detail Retention
figure
plot(newrealphase(:,300),'--')
hold on
scatter(1:512,IoutDCT(:,300),'.','green')

%% Evaluation Indicators
% RMES
RMES = sqrt(mean((IoutDCT(:)-realphase(:)).^2))
% PSNR
PSNR = 20*log10(255/sqrt(mean((IoutDCT(:)-realphase(:)).^2)))
% ENL
enl = ENL(IoutDCT)
% SSI
ssl2 = SSI(newrealphase,IoutDCT)
% NC
nc2= NC(realphase,IoutDCT)
% SSIM
sslm = SSIM(IoutDCT,realphase)
% FSIM
FSIM = FeatureSIM(realphase,IoutDCT)
