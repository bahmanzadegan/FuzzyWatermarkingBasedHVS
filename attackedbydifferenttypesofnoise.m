%% Image Watermarked attached by different types of noise.
[filename, pathname] = uigetfile({'*.tif','FIF file'},'Please select an image');
QueryPath=[pathname,filename];
Z=imread(QueryPath);

 S = imnoise(Z,'gaussian',0,0.001);
 D = imnoise(Z,'salt & pepper',0.001);
 A = imnoise(Z,'speckle',0.001);
 
% PSNR_1
[Mark_Extarct]=WatermarkDetection(S,Mark);
[SNR,PSNR] = PSNRfunction(Z,S);
PSNR_1 = PSNR;

% PSNR_2
[Mark_Extarct]=WatermarkDetection(D,Mark);
[SNR,PSNR] = PSNRfunction(Z,D);
PSNR_2 = PSNR;

% PSNR_3
[Mark_Extarct]=WatermarkDetection(A,Mark);
[SNR,PSNR] = PSNRfunction(Z,A);
PSNR_3 = PSNR;

% Plotting
figure,subplot(1,3,1),imshow(S),title({'Gaussian nois';['Var = ',num2str(0.001)];['PSNR = ',num2str(PSNR_1),' dB']});
subplot(1,3,2),imshow(D),title({'salt & pepper noise';['Var = ',num2str(0.001)];['PSNR = ',num2str(PSNR_2),' dB']});
subplot(1,3,3),imshow(A),title({'speckle noise';['Var = ',num2str(0.001)];['PSNR = ',num2str(PSNR_3),' dB']});

