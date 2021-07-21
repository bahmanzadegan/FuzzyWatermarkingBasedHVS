%function [Mark_Extarct]= WatermarkDetection(Watermarked,Mark)

%if (Watermarked == zeros )
    [filename, pathname] = uigetfile({'*.tif','FIF file'},'Please select Watemarked-image');
    QueryPath=[pathname,filename];
    Watermarked=imread(QueryPath);
    % imshow (QueryImage);
    disp('         Ok , Watermarked-image selected.');
    disp(' ')

    if ndims(Watermarked)==3
        Watermarked=rgb2gray(Watermarked);
    end
%end

V=blkproc(Watermarked,[8 8],@dct2); % enteghale tasvir be hozeye dct

u=1;
[x y]=size(Watermarked);
i=1;
j=1;

for i=i:8:x
    for j=j:8:y
        V_u = V(i:(i-1)+8,j:(j-1)+8); % mohasebeye har belak
        [t p]=size(V_u);
        for k=1:1:t
            for h=1:1:p
                V_1(u)=V_u(k,h); % chidane belakhaye poshte sare ham be soorate yek arraye 1bodi
                u=u+1;
            end
        end
    end
    j=1;
end

n=1;
load('D:\Uni901\Fuzzy\Project\Source code\fuzzy\Result\table');
[x y]=size(table);
for i=2:1:x
    for j=1:1:y-1
        B(1,j)=cast(table{i,j},'int32');
    end
    C1=V_1(B(1,j-2));
    C2=V_1(B(1,j-1));
    alfa = B(1,j);
    if C1>C2
        N(n)=0;
        n=n+1;
    else
        N(n)=1;
        n=n+1;
    end
end

[x y]=size(N);
i=1;
j=1;
m=0;
u=0;
h=1;
v=0;
l=1;
while h<=15 && l<=15
for i=1:1:x
    for j=j:1:15+m
        u=u+1;
        N_1(u)=N(i,j);
    end
    for h=h:1:15+v
           for l=1:1:15
           Mark_Extarct(h,l)=N_1(l); % chidane andis haye Qtable poshte sare ham be soorate yek arraye 1bodi
           end
           u=0;
           break;
    end
    m=m+15;
    j=j+1;
    h=h+1;
    v=v+15;
end
end

% Correlation value
Corelation_value = corr2(Mark,Mark_Extarct);
% Plotting
figure,imshow(Mark_Extarct),title({'Mark-Extarcted';['Correlation = ',num2str(Corelation_value)]});
imwrite(Mark_Extarct,'D:\Uni901\Fuzzy\Project\Source code\fuzzy\Result\Mark_Extarct.bmp','bmp');
