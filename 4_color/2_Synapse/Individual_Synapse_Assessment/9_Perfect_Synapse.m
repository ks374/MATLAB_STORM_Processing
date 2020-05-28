readpath = 'Z:\Chenghang\0_For_Processing\chenghaz_011_Sample_03_P4_2\analysis\elastic_align\Result_Syn\';
outpath = 'Z:\Chenghang\0_For_Processing\chenghaz_011_Sample_03_P4_2\analysis\elastic_align\Result_Syn\';
voxel = [15.5 15.5 70];
%
load([readpath 'G_paired_2.mat'])
load([readpath 'R_paired_2.mat'])
statsGwater = statsGwater_sss;
statsRwater = statsRwater_sss;
clear statsGwater_sss statsRwater_sss
%
expfolder = 'Z:\Chenghang\0_For_Processing\chenghaz_011_Sample_03_P4_2\analysis\elastic_align\storm_merged\';
files = [dir([expfolder '*.tif']) dir([expfolder '*.png'])];
info = imfinfo([expfolder files(1,1).name]);
num_images = numel(files);
furtherds = 1;

new_G = zeros(ceil((info(1,1).Height*furtherds)),ceil((info(1,1).Width*furtherds)),num_images,'uint8');
new_R = zeros(ceil((info(1,1).Height*furtherds)),ceil((info(1,1).Width*furtherds)),num_images,'uint8');
%%
statsG_temp = statsGwater;
for i = 1:size(statsG_temp,1)
    disp(i)
    PixelList = statsG_temp(i).PixelList;
    PixelValues = statsG_temp(i).PixelValues;
    for j = 1:size(PixelList,1)
        x = PixelList(j,2);
        y = PixelList(j,1);
        z = PixelList(j,3);
        new_G(x,y,z) = PixelValues(j);
    end
end
clear statsG_temp
statsG_temp = statsRwater;
for i = 1:size(statsG_temp,1)
    disp(i)
    PixelList = statsG_temp(i).PixelList;
    PixelValues = statsG_temp(i).PixelValues;
    for j = 1:size(PixelList,1)
        x = PixelList(j,2);
        y = PixelList(j,1);
        z = PixelList(j,3);
        new_R(x,y,z) = PixelValues(j);
    end
end
%%
for i=numel(statsGwater):-1:1

   minpix = min(statsGwater(i).PixelList);  maxpix = max(statsGwater(i).PixelList);
   min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
   max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end     
   BPp(i).mat = new_G(min2:max2,min1:max1,min3:max3); 
   BP2p(i).mat = new_R(min2:max2,min1:max1,min3:max3); 
end
%%
result = zeros(numel(statsGwater),1);
for jj=1:numel(statsGwater)
% jj = 1004;
   disp(['Calculating the angle of ' sprintf('%05d',jj)]);
   minpix = min(statsGwater(jj).PixelList);  maxpix = max(statsGwater(jj).PixelList);
   min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
   max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end        

   size1 = max1-min1 + 1; size2 = max2-min2 + 1; size3 = max3-min3 + 1;
   curr2 = false(size2, size1, size3);

   for j=1: numel(statsGwater(jj).PixelList(:,1))
       curr2(statsGwater(jj).PixelList(j,2)-min2+1, ...
          statsGwater(jj).PixelList(j,1)-min1+1, ...
          statsGwater(jj).PixelList(j,3)-min3+1)= 1;  
   end
%    curr1a = BPp(jj).mat;
   curr1b = BP2p(jj).mat;
   Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);

   size2= 105;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   
    Aa = permute(curr2,[1 3 2]);
    Ab = permute(curr1b,[1 3 2]);
    Ac = permute(curr4b,[1 3 2]);
clear A2a A2b A2c
    for k = 1:ceil(size(Aa,3))
     A2a(:,:,k) = imresize(Aa(:,:,k), [ceil(size(Aa,1)) ...
         ceil(size(Aa,2)*70/15.5)]);
     A2b(:,:,k) = imresize(Ab(:,:,k), [ceil(size(Ab,1)) ...
         ceil(size(Ab,2)*70/15.5)]);
     A2c(:,:,k) = imresize(Ac(:,:,k), [ceil(size(Ac,1)) ...
         ceil(size(Ac,2)*70/15.5)]);
    end
    curr2 = permute(A2a,[1 3 2]);
    curr1b = permute(A2b,[1 3 2]);
    curr4b = permute(A2c,[1 3 2]);
   
%       
   size1 = size(curr2,2); size2 = size(curr2,1); size3 = size(curr2,3);

   cg1.Connectivity = 26;
   cg1.ImageSize = [size2,size1,size3];
   cg1.NumObjects = 1;
   c = num2cell(find(curr2),1);
   cg1.PixelIdxList = c;
   
   cg2.Connectivity = 26;
   cg2.ImageSize = [size2,size1,size3];
   cg2.NumObjects = 1;
   temp_list = ones(numel(curr1b),1);
   temp_list(find(curr4b)) = 0;
   curr3 = curr1b;
   curr3(find(temp_list)) = 0;
   c = num2cell(find(curr3),1);
   cg2.PixelIdxList = c;
   
   ccg1 = regionprops3(cg1,'Orientation');
   ccg2 = regionprops3(cg2,'Orientation');

   angle1 = table2array(ccg1);
   angle2 = table2array(ccg2);
   angle1 = angle1./180.*3.14159;
   angle2 = angle2./180.*3.14159;
   result(jj) = cos(angle1(1))*cos(angle1(3))*cos(angle2(1))*cos(angle2(3)) + cos(angle1(3))*sin(angle1(1))*cos(angle2(3))*sin(angle2(1)) + sin(angle1(3)) * sin(angle2(3));

end
% result = abs(result);
%%
for i = 1:size(curr4b,3)
    imwrite(uint8(curr4b(:,:,i)),[outpath sprintf('%03d',i) '.tif']);
end
%%
for i = 1:size(curr2,3)
    imwrite(uint8(curr2(:,:,i)),[outpath sprintf('%03d',i) '.tif']);
end
%%
hist(result,80);
%%

val1w = result;
val2w = log(1+[statsGwater.Area]);
val2w = val2w';

Xn=80; Yn=80;
Xrange=[min(val1w) max(val1w)]; Yrange=[min(val2w) max(val2w)];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ; 

figure;
H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
close

cutoffg = 20;
H2 = H;
H2(H2>cutoffg)=cutoffg;

figure;
pcolor(X,Y,H2)
%%
% manually draw polygon on figure
currpoly=impoly
% return polygon coordinates
synapse_regiong=currpoly.getPosition

%Return centroid and stats lists for all clusters in selected polygon area
statsGwater_sss = statsGwater(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

statsGwater_ssn = statsGwater(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

size(statsGwater_sss,1)
%%
idx = result > 0.3;
% idx = result > 0.5;
statsGwater_sss = statsGwater(idx);
statsGwater_ssn = statsGwater(~idx);