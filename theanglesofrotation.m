%% Correlation value according to the angles of rotation

[filename, pathname] = uigetfile({'*.tif','FIF file'},'Please select Watemarked-image');
QueryPath=[pathname,filename];
Watemarked=imread(QueryPath);
% imshow (QueryImage);
disp('         Ok , Watermarked-image selected.');
disp(' ')

if ndims(Watemarked)==3
    Watemarked=rgb2gray(Watemarked);
end

%The Rotation Angle 1
Watemarked_Rotate = imrotate(Watemarked,1,'bicubic');
Watemarked_Rotate = imresize(Watemarked_Rotate,[256 256]);
[Mark_Extarct]=WatermarkDetection(Watemarked_Rotate,Mark);

%The Rotation Angle 2
Watemarked_Rotate = imrotate(Watemarked,2,'bicubic');
Watemarked_Rotate = imresize(Watemarked_Rotate,[256 256]);
[Mark_Extarct]=WatermarkDetection(Watemarked_Rotate,Mark);

%The Rotation Angle 3
Watemarked_Rotate = imrotate(Watemarked,3,'bicubic');
Watemarked_Rotate = imresize(Watemarked_Rotate,[256 256]);
[Mark_Extarct]=WatermarkDetection(Watemarked_Rotate,Mark);

%The Rotation Angle 4
Watemarked_Rotate = imrotate(Watemarked,4,'bicubic');
Watemarked_Rotate = imresize(Watemarked_Rotate,[256 256]);
[Mark_Extarct]=WatermarkDetection(Watemarked_Rotate,Mark);

%The Rotation Angle 5
Watemarked_Rotate = imrotate(Watemarked,5,'bicubic');
Watemarked_Rotate = imresize(Watemarked_Rotate,[256 256]);
[Mark_Extarct]=WatermarkDetection(Watemarked_Rotate,Mark);
