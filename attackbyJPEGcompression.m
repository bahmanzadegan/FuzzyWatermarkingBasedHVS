
[filename, pathname] = uigetfile({'*.tif','FIF file'},'Please select Watemarked-image');
QueryPath=[pathname,filename];
Watemarked=imread(QueryPath);
% imshow (QueryImage);
disp('         Ok , Watermarked-image selected.');
disp(' ')

if ndims(Watemarked)==3
    Watemarked=rgb2gray(Watemarked);
end
I = Watemarked ;

% JPEG compression
jpegcompression(Watemarked ,'lena_compressed.mat');
img = Watemarked;
% JPEG Watermarking
[Mark_Extarct]= WatermarkDetection(Watemarked,Mark);
img_1 = Watemarked ;
% System performances
[CR,BPP,PSNR,MSE,SNR] = systemperformances(img,I,'lena_compressed.mat');

% Plotting
figure,subplot(1,2,1),imshow(I),title({'Watemarked image';'';''});
subplot(1,2,2),imshow(img),title({'JPEG Compressed image';['PSNR = ',num2str(PSNR),' dB'];['Compression ratio = ',num2str(CR)]});