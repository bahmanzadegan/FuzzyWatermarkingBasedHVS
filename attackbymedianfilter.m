%% Correlation values after attack by median filter.

[filename, pathname] = uigetfile({'*.tif','FIF file'},'Please select Watemarked-image');
QueryPath=[pathname,filename];
Watemarked=imread(QueryPath);
% imshow (QueryImage);
disp('         Ok , Watermarked-image selected.');
disp(' ')

if ndims(Watemarked)==3
    Watemarked=rgb2gray(Watemarked);
end


%Median filter size 3x3
Watemarked_filter  = medfilt2(Watemarked, [3 3]);
[Mark_Extarct]= WatermarkDetection(Watemarked_filter,Mark);

%Median filter size 5x5
Watemarked_filter  = medfilt2(Watemarked, [5 5]);
[Mark_Extarct]=WatermarkDetection(Watemarked_filter,Mark);

%Median filter size 7x7
Watemarked_filter  = medfilt2(Watemarked, [7 7]);
[Mark_Extarct]=WatermarkDetection(Watemarked_filter,Mark);