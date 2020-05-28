clear all;clc
%
base_path = 'Z:\Chenghang\chenghaz_015_B2_P8_1eye_one_area\';
exp_folder = [base_path 'analysis\'];
path = [exp_folder  'elastic_align\'];
mergedpath = [path 'storm_merged/'];
%
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);
disp(num_images)
%set voxel size in nm
voxel=[15.5, 15.5, 70];

inpath = [exp_folder 'Result\5_V_Syn\'];
outpath = [exp_folder 'Result\'];
%
load([inpath 'R_paired_3.mat']);
statsRwater = statsRwater_ssss;
%
load([inpath 'G_paired_3.mat']);
statsGwater = statsGwater_ssss;

clear statsRwater_ssss statsGwater_ssss
%
clear BP BG bg2
disp('allocating arrays')
BP2 = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsGwater)
    disp(int2str(i))
    BP2(statsGwater(i).PixelIdxList)=statsGwater(i).PixelValues;  
end
%
for i=numel(statsRwater):-1:1
   minpix = min(statsRwater(i).PixelList);  maxpix = max(statsRwater(i).PixelList);
   min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
   max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end      
   BP2p(i).mat = BP2(min2:max2,min1:max1,min3:max3); 
end
%
%If skip this step, most the closest distances will be 0. So there a
%threshold is used on the statsGwater (BP2), so that only the core region
%of the cluster will be calculated. 
for i = 1:1:numel(BP2p)
    temp = BP2p(i).mat(find(BP2p(i).mat));
    thre = multithresh(temp,2);
    thre = thre(2);
    BP2p(i).mat(BP2p(i).mat <= thre) = 0;
end
%
Dist = zeros(numel(statsRwater),1);
parfor jj=1:numel(statsRwater)
   disp(jj)
   minpix = min(statsRwater(jj).PixelList);  maxpix = max(statsRwater(jj).PixelList);
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
   for j=1: numel(statsRwater(jj).PixelList(:,1))
       curr2(statsRwater(jj).PixelList(j,2)-min2+1, ...
          statsRwater(jj).PixelList(j,1)-min1+1, ...
          statsRwater(jj).PixelList(j,3)-min3+1)= 1;  
   end
   curr1b = BP2p(jj).mat;
   Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
   Dist(jj) = min(Dg(find(curr1b)));
   %Dg_temp used to verify the code works find. It works fine indeed... 
   %Dg_temp(find(curr1b)) = Dg(find(curr1b));
end
%
if ~exist([outpath '6_Dist\'])
    mkdir([outpath '6_Dist\'])
end
save([outpath '6_Dist\Dist_edge.mat'],'Dist');