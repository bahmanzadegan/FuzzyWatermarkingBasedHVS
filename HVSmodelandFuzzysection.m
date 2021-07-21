function [alfa,Q_1] = HVSmodelandFuzzysection(V_dc,V_1,k)
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
mfedit(F)
%plotfis(F)
Fuzzy_memberships = evalfis([L ; T_k_1],F); %defuzz (Centroid)

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