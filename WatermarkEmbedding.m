clc;
clear all;

[filename, pathname] = uigetfile({'*.tif','FIF file'},'Please select Orginal-image');
QueryPath=[pathname,filename];
QueryImage=imread(QueryPath);
% imshow (QueryImage);
disp('         Ok , Orginal-image selected.');
disp(' ')

if ndims(QueryImage)==3
    QueryImage=rgb2gray(QueryImage);
end

I = imresize(QueryImage,[256 256]);
imwrite(I,'D:\Uni901\Fuzzy\Project\Source code\fuzzy\Result\Orginal_image.tif','tif');
%figure,imshow(I); % tasvire vorudi

V=blkproc(I,[8 8],@dct2); % enteghale tasvir be hozeye dct

u=1;
[x y]=size(I);
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

[x y]=size(I);
m=x;
b=y;
m = m/8;
b = b/8;
xy=m*b; % mohasebeye tedade black'ha

k=1; % tedade belak'ha

    for i=1:8:x 
        for j=1:8:y
           V_dc(k)=V(i,j); % zazribe DC black k'om
           k=k+1;
        end
    end
%% function (HVSmodelandFuzzysection)
%function [alfa,Q_1] = HVSmodelandFuzzysection(V_dc,V_1,k)
%% L(k)
 Vmean=0; 
    for t=1:1:k-1
        Vmean = Vmean + V_dc(t); % mohasebeye miyangine zarayebe DC
    end
    Vmean = Vmean / k-1 ;    
 
 r = 0.649; % khode maghale pishnahad dade ast
 
    for f=1:1:k-1
        L(f)= (V_dc(f) / Vmean) .^ r; %Luminance sensitivity : L(k)
    end
%% Quantization table

Qt = [16	11	10	16	24	40	51	61;
      12	12	14	19	26	58	60	55;
      14	13	16	24	40	57	69	56;
      14	17	22	29	51	87	80	62;
      18	22	37	56	68	109	103	77;
      24	35	55	64	81	104	113	92;
      49	64	78	87	103	121	120	101;
      72	92	95	98	112	100	103	99 ]; % JPEG quantization table (luminance)

u=1;
[x y]=size(Qt);
i=1;
j=1;

for i=i:1:x
    for j=j:1:y
        Q_u = Qt(i:(i-1)+1,j:(j-1)+1); % mohasebeye har andise quantization table
        [t p]=size(Q_u);
        for k=1:1:t
            for h=1:1:p
                Q_1(u)=Q_u(k,h); % chidane andis haye Qtable poshte sare ham be soorate yek arraye 1bodi
                u=u+1;
            end
        end
    end
    j=1;
end
%% T(k)
[e b]=size(V_1);
[d g]=size(Q_1);
f=0;
x=1;
T_k = 0;
 
while f < b
for i=1:1:e
    for j=1:1:g
        if f==0
        T = cond(round(V_1(i,j) / Q_1(i,j))) ; % mohasebeye round kardane adad va dar sorate gheyre 0 shodan,meghdare 1 gereftan T
            if T == Inf
               T=0;
            end
        else
        T = cond(round(V_1(i,j+(63*x)) / Q_1(i,j))) ;
            if T == Inf
               T=0;
            end
        end
        T_k = T_k + T ; % mohasebeye majmooye andishaye gheyre 0 dar in block
    end
    f=f+64;
end
T_k_1(x)=T_k; % mohasebeye T(k)
x = x+1 ;
T_k = 0 ;
end

%% Fis
F = readfis('D:\Uni901\Fuzzy\Project\Source code\fuzzy\Fuzzy(File.fis)\fis1.fis'); 
%mfedit(F)
%plotfis(F)
Fuzzy_memberships = evalfis([T_k_1 ; L],F); %defuzz (Centroid)

%% mohasebeye alfa
[x y]=size(Fuzzy_memberships);
for i=1:1:x
    for j=1:1:y
        Fuzzy_memberships_1(j,i) = Fuzzy_memberships(i,j); % tabdile array (:,1) be arraye (1,:)
    end
end
[z e]=size(Fuzzy_memberships_1);
[n m]=size(Q_1);

w=1;
j=1;
while w<=e
for i=1:1:n
    for j=j:1:m
        alfa(i,w) = Fuzzy_memberships_1(1,w) * Q_1(i,j); % zarbe khoro0jie defuzzy ba frequency sensivitity
    end
    k=k+64;
    w=w+1;
end
end
%[alfa,Q_1] = HVSmodelandFuzzysection (V_dc,V_1,k);

%% The Watermark Embededding Process
%disp('Step_2 : ( Selecting Mark )')
[filename, pathname] = uigetfile({'*.bmp','Bit Map file'},'Please select a Mark');
MarkPath=[pathname,filename];
Mark=imread(MarkPath);

disp('         Ok , Mark selected.');
disp(' ')

if ndims(Mark)==3
    Mark=rgb2gray(Mark);
end

Mark = imresize(Mark,[15 15]);
imwrite(Mark,'D:\Uni901\Fuzzy\Project\Source code\fuzzy\Result\Orginal_mark.bmp','bmp');
% figure,
% imshow(Mark);

u=0;
[x y]=size(Mark);
i=1;
j=1;
for i=1:1:x
    for j=j:1:y
        u=u+1;
       Mark_1(u)=Mark(i,j); % chidane andis haye Qtable poshte sare ham be soorate yek arraye 1bodi
    end
    j=1;
end

for i=1:1:u
    if Mark_1(i)> 0
       Mark_1(i)=1; % chidane andis haye Qtable poshte sare ham be soorate yek arraye 1bodi
    else 
       Mark_1(i)=0;
    end
end

%Mark_1=[0,1,1,0,1,1,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0]; % mark
  [r p]=size(Mark_1);
  
  h=1;
  l=1;
  o=1;
  [a k]=size(Q_1);
  table = {'C1-index','C2-index','Block','Quantization'};
  w=2;
  [x y]=size(V_1);
  f=1;
  t=1;
  
  while o<=p
  Qt_index_1=0;
  Qt_index=-1;
  C1=0;
  C2=0;
  d=0;
    if Mark_1(o)==1 % mohasebeye 2 andise ba maghadire yeksan dar jadvale quantization table
      
      [x y]=size(V_1); % mohasebeye C1,C2
      for u=1:1:x
      for f=f:64:y
           
         for i=1:1:a
             for j=2:1:k
                Qt_index=Q_1(i,j);
                
                for z=1:1:a
                     for s=2:1:k
                         if i==z && s==j
                         else
                              Qt_index_1=Q_1(z,s);
                              if Qt_index_1==Qt_index
                                  
                                  
                                          C1=V_1(1,(f-1)+j);
                                          C2=V_1(1,(f-1)+s);
                                          if (f==1)
                                              alfa_1 = alfa(u,1);
                                          else
                                              alfa_1 = alfa(u,(nearest((f-1)/64))+1);
                                          end
                                          d=0;
                                          if (C1<C2 && C1-C2<=alfa_1) % niaz be chek shodane (C1~=C1's_dataset && C2~=C2's_dataset)
                                                 if o ==1
                                                    table(w,:) = {(f-1)+j,(f-1)+s,(nearest((f/64)+0.5)),Qt_index_1};
                                                    w=w+1;
                                                    C1_1 = C1 - (alfa_1 / 2);
                                                    C2_2 = C2 + (alfa_1 / 2);
                                                    
                                                    V_1(1,(f-1)+j) = C1_1;
                                                    V_1(1,(f-1)+s) = C2_2;
                                                    
                                                 else
                                                    for g=2:1:o
                                                       if  (table{g,1} == (f-1)+j && table{g,2} == (f-1)+s  && table{g,4} == Qt_index_1) || (table{g,3} == nearest((f/64)+0.5)) 
                                                           d=d+1;
                                                       end
                                                    end
                                                    if d~=0
                                                       %break;
                                                    else
                                                       table(w,:) = {(f-1)+j,(f-1)+s,(nearest((f/64)+0.5)),Qt_index_1};
                                                       w=w+1;
                                                       
                                                       C1_1 = C1 - (alfa_1 / 2);
                                                       C2_2 = C2 + (alfa_1 / 2);
                                                    
                                                       V_1(1,(f-1)+j) = C1_1;
                                                       V_1(1,(f-1)+s) = C2_2;
                                                   %darje market (formule(5))
                                                   %dataset = V(h*i,l*j), V(h*z,l*s) ,Qt_index , o (jahate estefade az an dar estekhraje market) 
                                                    end
                                                 end
                                           end
                              if Qt_index_1==Qt_index && (C1<C2 && C1-C2<=alfa_1) && d==0 % niaz be chek shodane (C1~=C1's_dataset && C2~=C2's_dataset)
                                    break;
                              end
                              end 
                         end
                     if Qt_index_1==Qt_index && (C1<C2 && C1-C2<=alfa_1) && d==0 % niaz be chek shodane (C1~=C1's_dataset && C2~=C2's_dataset)
                           break;
                     end
                     end
                end
              if Qt_index_1==Qt_index && (C1<C2 && C1-C2<=alfa_1)
                     break;
              end
              end
          if Qt_index_1==Qt_index && (C1<C2 && C1-C2<=alfa_1)
                break;
          end
          end
       if Qt_index_1==Qt_index && (C1<C2 && C1-C2<=alfa_1)
             break;
       end
       end
       if Qt_index_1==Qt_index && (C1<C2 && C1-C2<=alfa_1)
             break;
       end
       end     
         
    elseif Mark_1(o)==0
       [x y]=size(V_1); % mohasebeye C1,C2
       for u=1:1:x
       for t=t:64:y
            
         for i=1:1:a
             for j=2:1:k
                Qt_index=Q_1(i,j);
                
                for z=1:1:a
                     for s=2:1:k
                         if i==z && s==j
                         else
                              Qt_index_1=Q_1(z,s);
                              if Qt_index_1==Qt_index
                                  
                                  
                                          C1=V_1(1,(t-1)+j);
                                          C2=V_1(1,(t-1)+s);
                                          if (t==1)
                                              alfa_1 = alfa(u,1);
                                          else
                                              alfa_1 = alfa(u,(nearest((t-1)/64))+1);
                                          end
                                          d=0;
                                          if (C1>C2 && C1-C2<=alfa_1 ) % niaz be chek shodane (C1~=C1's_dataset && C2~=C2's_dataset)
                                                 if o ==1
                                                    table(w,:) = {(t-1)+j,(t-1)+s,(nearest((t/64)+0.5)),Qt_index_1};
                                                    w=w+1;
                                                    C1_1 = C1 + (alfa_1 / 2);
                                                    C2_2 = C2 - (alfa_1 / 2);
                                                    
                                                    V_1(1,(t-1)+j) = C1_1;
                                                    V_1(1,(t-1)+s) = C2_2;
                                                 else
                                                    for g=2:1:o
                                                       if  (table{g,1} == (t-1)+j && table{g,2} == (t-1)+s && table{g,4} == Qt_index_1) || (table{g,3} == nearest((t/64)+0.5)) 
                                                           d=d+1;
                                                       end
                                                    end
                                                    if d~=0
                                                       %break;
                                                    else
                                                       table(w,:) = {(t-1)+j,(t-1)+s,(nearest((t/64)+0.5)),Qt_index_1};
                                                       w=w+1;
                                                       C1_1 = C1 + (alfa_1 / 2);
                                                       C2_2 = C2 - (alfa_1 / 2);
                                                    
                                                       V_1(1,(t-1)+j) = C1_1;
                                                       V_1(1,(t-1)+s) = C2_2;
                                                   %darje market (formule(5))
                                                   %dataset = V(h*i,l*j), V(h*z,l*s) ,Qt_index , o (jahate estefade az an dar estekhraje market) 
                                                    end
                                                 end
                                           end
                              if Qt_index_1==Qt_index && (C1>C2 && C1-C2<=alfa_1 )&& d==0 % niaz be chek shodane (C1~=C1's_dataset && C2~=C2's_dataset)
                                     break;
                              end 
                              end 
                         end
                     if Qt_index_1==Qt_index && (C1>C2 && C1-C2<=alfa_1 ) && d==0 % niaz be chek shodane (C1~=C1's_dataset && C2~=C2's_dataset)
                          break;
                     end 
                     end
                end
              if Qt_index_1==Qt_index && (C1>C2 && C1-C2<=alfa_1 )
                   break;
              end
              end
          if Qt_index_1==Qt_index && (C1>C2 && C1-C2<=alfa_1 )
                break;
          end
          end
        if Qt_index_1==Qt_index && (C1>C2 && C1-C2<=alfa_1 )
              break;
        end
        end
        if Qt_index_1==Qt_index && (C1>C2 && C1-C2<=alfa_1 )
              break;
        end
        end
 
    end %end if
      
    if (((C1<C2 && C1-C2<=alfa_1) && Mark_1(o)==1 ) || ((C1>C2 && C1-C2<=alfa_1 )&& Mark_1(o)==0))
      o=o+1;
    end

  end

u=0;
l=1;
h=1;
n=8;
m=8;
d=0;
b=1;
while (n <= 256 && m <= 256)
        for k=l:1:n
            for h=h:1:m
                u=u+1;
                V_u(k,h)=V_1(u); % chidane belakhaye poshte sare ham
            end
            h=h-7;
            
        end
        d=d+1;
        if (d==32)
         l=l+8;
         n=n+8;
        h=1;
        if u == 65536
            m=257;
        else
        m=8;
        end
        d=0;
        b=b+1;
        else
        h=h+8;
        m=m+8;
        end
end
Z=blkproc(V_u,[8 8],@idct2);
Z=uint8(Z);
%figure,imshow(Z);
save('D:\Uni901\Fuzzy\Project\Source code\fuzzy\Result\table','table'); % eror darad
imwrite(Z,'D:\Uni901\Fuzzy\Project\Source code\fuzzy\Result\Watermarked_image.tif','tif');

%% Performance evaluation

% PSNR Function
[SNR,PSNR] = PSNRfunction(I,Z);
disp('PSNR Peak Signal to Noise Ratio');
disp(PSNR);

% Correlation value
Corelation_value = corr2(I,Z);
disp('Correlation value');
disp(Corelation_value);

% Bit error rate
% disp('Bit error rate');
% disp(nnz(message-rec_M)/numel(M));

%% Plotting
figure,subplot(1,2,1),imshow(I),title({'Orginal medical image';'';''});
subplot(1,2,2),imshow(Z),title({'Medical image watermarked';['PSNR = ',num2str(PSNR),' dB'];['Correlation = ',num2str(Corelation_value)]});
