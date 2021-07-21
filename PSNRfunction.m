
function [SNR,PSNR] = PSNRfunction(OrginalImage,WatermarkedImage)

I0     = double(OrginalImage);

I1     = double(WatermarkedImage);

Id     = (I0-I1);

signal = sum(sum(I0.^2));

noise  = sum(sum(Id.^2));

SNR    = 10*log10(signal/noise);

disp(SNR);

MSE  = noise/numel(I0);

peak = max(I0(:));

PSNR = 10*log10(peak^2/MSE);