%%
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

outpath = [exp_folder 'Result\'];
outpath_Vsyn = [outpath '5_V_Syn\'];
%
load([outpath_Vsyn 'R_paired_3.mat']);
statsRwater = statsRwater_ssss;
clear statsGwater
%
clear statsGa2s statsGwater
load([outpath_Vsyn 'G_paired_3.mat']);
statsGwater = statsGwater_ssss;
%
% num_images = num_images - 2;
clear BP BG bg2
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
BP2 = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsGwater)
    disp(int2str(i))
    BP(statsGwater(i).PixelIdxList)=statsGwater(i).PixelValues;  
end
%
%BP2 = BP;
for i = 1:numel(statsRwater)
    disp(int2str(i))
    BP2(statsRwater(i).PixelIdxList)=statsRwater(i).PixelValues;  
end
%
disp('loading data')
%
%if not present, load in stats list

%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsGwater):-1:1
    disp(i)
    
statsGwater(i).tints_p5 = [];
statsGwater(i).volume_p5 = [];
statsGwater(i).area_p5 = [];
statsGwater(i).WeightedCentroid_p5 = [];
centG_p5(i,1:3) = 0;
tintsG_p5(i) = 0;
volumeG_p5(i) =0;

statsGwater(i).tints_p10 = [];
statsGwater(i).volume_p10 = [];
statsGwater(i).area_p10 = [];
statsGwater(i).WeightedCentroid_p10 = [];
centG_p10(i,1:3) = 0;
tintsG_p10(i) = 0;
volumeG_p10(i) =0;

statsGwater(i).tints_p15 = [];
statsGwater(i).volume_p15 = [];
statsGwater(i).area_p15 = [];
statsGwater(i).WeightedCentroid_p15 = [];
centG_p15(i,1:3) = 0;
tintsG_p15(i) = 0;
volumeG_p15(i) =0;

statsGwater(i).tints_p20 = [];
statsGwater(i).volume_p20 = [];
statsGwater(i).area_p20 = [];
statsGwater(i).WeightedCentroid_p20 = [];
centG_p20(i,1:3) = 0;
tintsG_p20(i) = 0;
volumeG_p20(i) =0;

statsGwater(i).tints_p25 = [];
statsGwater(i).volume_p25 = [];
statsGwater(i).area_p25 = [];
statsGwater(i).WeightedCentroid_p25 = [];
centG_p25(i,1:3) = 0;
tintsG_p25(i) = 0;
volumeG_p25(i) =0;


statsGwater(i).tints_p30 = [];
statsGwater(i).volume_p30 = [];
statsGwater(i).area_p30 = [];
statsGwater(i).WeightedCentroid_p30 = [];
centG_p30(i,1:3) = 0;
tintsG_p30(i) = 0;
volumeG_p30(i) =0;

statsGwater(i).tints_p35 = [];
statsGwater(i).volume_p35 = [];
statsGwater(i).area_p35 = [];
statsGwater(i).WeightedCentroid_p35 = [];
centG_p35(i,1:3) = 0;
tintsG_p35(i) = 0;
volumeG_p35(i) =0;

statsGwater(i).tints_p40 = [];
statsGwater(i).volume_p40 = [];
statsGwater(i).area_p40 = [];
statsGwater(i).WeightedCentroid_p40 = [];
centG_p40(i,1:3) = 0;
tintsG_p40(i) = 0;
volumeG_p40(i) =0;

statsGwater(i).tints_p45 = [];
statsGwater(i).volume_p45 = [];
statsGwater(i).area_p45 = [];
statsGwater(i).WeightedCentroid_p45 = [];
centG_p45(i,1:3) = 0;
tintsG_p45(i) = 0;
volumeG_p45(i) =0;

statsGwater(i).tints_p50 = [];
statsGwater(i).volume_p50 = [];
statsGwater(i).area_p50 = [];
statsGwater(i).WeightedCentroid_p50 = [];
centG_p50(i,1:3) = 0;
tintsG_p50(i) = 0;
volumeG_p50(i) =0;

statsGwater(i).tints_p55 = [];
statsGwater(i).volume_p55 = [];
statsGwater(i).area_p55 = [];
statsGwater(i).WeightedCentroid_p55 = [];
centG_p55(i,1:3) = 0;
tintsG_p55(i) = 0;
volumeG_p55(i) =0;

statsGwater(i).tints_p60 = [];
statsGwater(i).volume_p60 = [];
statsGwater(i).area_p60 = [];
statsGwater(i).WeightedCentroid_p60 = [];
centG_p60(i,1:3) = 0;
tintsG_p60(i) = 0;
volumeG_p60(i) =0;

statsGwater(i).tints_p65 = [];
statsGwater(i).volume_p65 = [];
statsGwater(i).area_p65 = [];
statsGwater(i).WeightedCentroid_p65 = [];
centG_p65(i,1:3) = 0;
tintsG_p65(i) = 0;
volumeG_p65(i) =0;

statsGwater(i).tints_p70 = [];
statsGwater(i).volume_p70 = [];
statsGwater(i).area_p70 = [];
statsGwater(i).WeightedCentroid_p70 = [];
centG_p70(i,1:3) = 0;
tintsG_p70(i) = 0;
volumeG_p70(i) =0;
%try preproccessing the BP varialbe in a for loop
   minpix = min(statsGwater(i).PixelList);  maxpix = max(statsGwater(i).PixelList);
   min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
   max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end      
   BPp(i).mat = BP(min2:max2,min1:max1,min3:max3); 
   BP2p(i).mat = BP2(min2:max2,min1:max1,min3:max3); 
end
%

for jj=1:numel(statsGwater)
   %
   disp(jj)
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
   curr1a = BPp(jj).mat;
   curr1b = BP2p(jj).mat;
   Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
%
    size2= 5;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m5 = [x,y,z];
   RPixelValues_m5 = zeros(numel(PixelList_m5(:,1)),1);
   for kk=1:numel(PixelList_m5(:,1))
      RPixelValues_m5(kk) = curr1b(PixelList_m5(kk,2),...
          PixelList_m5(kk,1), PixelList_m5(kk,3));
   end

   statsGwater(jj).tints_p5 = sum([RPixelValues_m5]);
   statsGwater(jj).area_p5 = numel([RPixelValues_m5]);
   statsGwater(jj).volume_p5 = numel([RPixelValues_m5]>0);
   statsGwater(jj).WeightedCentroid_p5(1) = sum([PixelList_m5(:,1)].*...
           double([RPixelValues_m5]))/(sum([RPixelValues_m5]));
   statsGwater(jj).WeightedCentroid_p5(2) = sum([PixelList_m5(:,2)].*...
           double([RPixelValues_m5]))/(sum([RPixelValues_m5]));
   statsGwater(jj).WeightedCentroid_p5(3) =  sum([PixelList_m5(:,3)].*...
           double([RPixelValues_m5]))/(sum([RPixelValues_m5]));
%also add these statxture, 
   centG_p5(jj,:) = statsGwater(jj).WeightedCentroid_p5;
   tintsG_p5(jj) = statsGwater(jj).tints_p5;
   volumeG_p5(jj) =statsGwater(jj).volume_p5;
   
   
   size2= 10;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m10 = [x,y,z];
   RPixelValues_m10 = zeros(numel(PixelList_m10(:,1)),1);
   for kk=1:numel(PixelList_m10(:,1))
      RPixelValues_m10(kk) = curr1b(PixelList_m10(kk,2),...
          PixelList_m10(kk,1), PixelList_m10(kk,3));
   end

   statsGwater(jj).tints_p10 = sum([RPixelValues_m10]);
   statsGwater(jj).area_p10 = numel([RPixelValues_m10]);
   statsGwater(jj).volume_p10 = numel([RPixelValues_m10]>0);
   statsGwater(jj).WeightedCentroid_p10(1) = sum([PixelList_m10(:,1)].*...
           double([RPixelValues_m10]))/(sum([RPixelValues_m10]));
   statsGwater(jj).WeightedCentroid_p10(2) = sum([PixelList_m10(:,2)].*...
           double([RPixelValues_m10]))/(sum([RPixelValues_m10]));
   statsGwater(jj).WeightedCentroid_p10(3) =  sum([PixelList_m10(:,3)].*...
           double([RPixelValues_m10]))/(sum([RPixelValues_m10]));
%also add these statxture, 
   centG_p10(jj,:) = statsGwater(jj).WeightedCentroid_p10;
   tintsG_p10(jj) = statsGwater(jj).tints_p10;
   volumeG_p10(jj) =statsGwater(jj).volume_p10;

   size2= 15;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m15 = [x,y,z];
   RPixelValues_m15 = zeros(numel(PixelList_m15(:,1)),1);
   for kk=1:numel(PixelList_m15(:,1))
      RPixelValues_m15(kk) = curr1b(PixelList_m15(kk,2),...
          PixelList_m15(kk,1), PixelList_m15(kk,3));
   end

   statsGwater(jj).tints_p15 = sum([RPixelValues_m15]);
   statsGwater(jj).area_p15 = numel([RPixelValues_m15]);
   statsGwater(jj).volume_p15 = numel([RPixelValues_m15]>0);
   statsGwater(jj).WeightedCentroid_p15(1) = sum([PixelList_m15(:,1)].*...
           double([RPixelValues_m15]))/(sum([RPixelValues_m15]));
   statsGwater(jj).WeightedCentroid_p15(2) = sum([PixelList_m15(:,2)].*...
           double([RPixelValues_m15]))/(sum([RPixelValues_m15]));
   statsGwater(jj).WeightedCentroid_p15(3) =  sum([PixelList_m15(:,3)].*...
           double([RPixelValues_m15]))/(sum([RPixelValues_m15]));
%also add these statxture, 
   centG_p15(jj,:) = statsGwater(jj).WeightedCentroid_p15;
   tintsG_p15(jj) = statsGwater(jj).tints_p15;
   volumeG_p15(jj) =statsGwater(jj).volume_p15;
   
   size2= 20;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m20 = [x,y,z];
   RPixelValues_m20 = zeros(numel(PixelList_m20(:,1)),1);
   for kk=1:numel(PixelList_m20(:,1))
      RPixelValues_m20(kk) = curr1b(PixelList_m20(kk,2),...
          PixelList_m20(kk,1), PixelList_m20(kk,3));
   end

   statsGwater(jj).tints_p20 = sum([RPixelValues_m20]);
   statsGwater(jj).area_p20 = numel([RPixelValues_m20]);
   statsGwater(jj).volume_p20 = numel([RPixelValues_m20]>0);
   statsGwater(jj).WeightedCentroid_p20(1) = sum([PixelList_m20(:,1)].*...
           double([RPixelValues_m20]))/(sum([RPixelValues_m20]));
   statsGwater(jj).WeightedCentroid_p20(2) = sum([PixelList_m20(:,2)].*...
           double([RPixelValues_m20]))/(sum([RPixelValues_m20]));
   statsGwater(jj).WeightedCentroid_p20(3) =  sum([PixelList_m20(:,3)].*...
           double([RPixelValues_m20]))/(sum([RPixelValues_m20]));
%also add these statxture, 
   centG_p20(jj,:) = statsGwater(jj).WeightedCentroid_p20;
   tintsG_p20(jj) = statsGwater(jj).tints_p20;
   volumeG_p20(jj) =statsGwater(jj).volume_p20;
   
   size2= 25;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m25 = [x,y,z];
   RPixelValues_m25 = zeros(numel(PixelList_m25(:,1)),1);
   for kk=1:numel(PixelList_m25(:,1))
      RPixelValues_m25(kk) = curr1b(PixelList_m25(kk,2),...
          PixelList_m25(kk,1), PixelList_m25(kk,3));
   end

   statsGwater(jj).tints_p25 = sum([RPixelValues_m25]);
   statsGwater(jj).area_p25 = numel([RPixelValues_m25]);
   statsGwater(jj).volume_p25 = numel([RPixelValues_m25]>0);
   statsGwater(jj).WeightedCentroid_p25(1) = sum([PixelList_m25(:,1)].*...
           double([RPixelValues_m25]))/(sum([RPixelValues_m25]));
   statsGwater(jj).WeightedCentroid_p25(2) = sum([PixelList_m25(:,2)].*...
           double([RPixelValues_m25]))/(sum([RPixelValues_m25]));
   statsGwater(jj).WeightedCentroid_p25(3) =  sum([PixelList_m25(:,3)].*...
           double([RPixelValues_m25]))/(sum([RPixelValues_m25]));
%also add these statxture, 
   centG_p25(jj,:) = statsGwater(jj).WeightedCentroid_p25;
   tintsG_p25(jj) = statsGwater(jj).tints_p25;
   volumeG_p25(jj) =statsGwater(jj).volume_p25;
   
   
   size2= 30;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m30 = [x,y,z];
   RPixelValues_m30 = zeros(numel(PixelList_m30(:,1)),1);
   for kk=1:numel(PixelList_m30(:,1))
      RPixelValues_m30(kk) = curr1b(PixelList_m30(kk,2),...
          PixelList_m30(kk,1), PixelList_m30(kk,3));
   end

   statsGwater(jj).tints_p30 = sum([RPixelValues_m30]);
   statsGwater(jj).area_p30 = numel([RPixelValues_m30]);
   statsGwater(jj).volume_p30 = numel([RPixelValues_m30]>0);
   statsGwater(jj).WeightedCentroid_p30(1) = sum([PixelList_m30(:,1)].*...
           double([RPixelValues_m30]))/(sum([RPixelValues_m30]));
   statsGwater(jj).WeightedCentroid_p30(2) = sum([PixelList_m30(:,2)].*...
           double([RPixelValues_m30]))/(sum([RPixelValues_m30]));
   statsGwater(jj).WeightedCentroid_p30(3) =  sum([PixelList_m30(:,3)].*...
           double([RPixelValues_m30]))/(sum([RPixelValues_m30]));
%also add these statxture, 
   centG_p30(jj,:) = statsGwater(jj).WeightedCentroid_p30;
   tintsG_p30(jj) = statsGwater(jj).tints_p30;
   volumeG_p30(jj) =statsGwater(jj).volume_p30;
   
   size2= 35;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m35 = [x,y,z];
   RPixelValues_m35 = zeros(numel(PixelList_m35(:,1)),1);
   for kk=1:numel(PixelList_m35(:,1))
      RPixelValues_m35(kk) = curr1b(PixelList_m35(kk,2),...
          PixelList_m35(kk,1), PixelList_m35(kk,3));
   end

   statsGwater(jj).tints_p35 = sum([RPixelValues_m35]);
   statsGwater(jj).area_p35 = numel([RPixelValues_m35]);
   statsGwater(jj).volume_p35 = numel([RPixelValues_m35]>0);
   statsGwater(jj).WeightedCentroid_p35(1) = sum([PixelList_m35(:,1)].*...
           double([RPixelValues_m35]))/(sum([RPixelValues_m35]));
   statsGwater(jj).WeightedCentroid_p35(2) = sum([PixelList_m35(:,2)].*...
           double([RPixelValues_m35]))/(sum([RPixelValues_m35]));
   statsGwater(jj).WeightedCentroid_p35(3) =  sum([PixelList_m35(:,3)].*...
           double([RPixelValues_m35]))/(sum([RPixelValues_m35]));
%also add these statxture, 
   centG_p35(jj,:) = statsGwater(jj).WeightedCentroid_p35;
   tintsG_p35(jj) = statsGwater(jj).tints_p35;
   volumeG_p35(jj) =statsGwater(jj).volume_p35;
   
   
   size2= 40;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m40 = [x,y,z];
   RPixelValues_m40 = zeros(numel(PixelList_m40(:,1)),1);
   for kk=1:numel(PixelList_m40(:,1))
      RPixelValues_m40(kk) = curr1b(PixelList_m40(kk,2),...
          PixelList_m40(kk,1), PixelList_m40(kk,3));
   end

   statsGwater(jj).tints_p40 = sum([RPixelValues_m40]);
   statsGwater(jj).area_p40 = numel([RPixelValues_m40]);
   statsGwater(jj).volume_p40 = numel([RPixelValues_m40]>0);
   statsGwater(jj).WeightedCentroid_p40(1) = sum([PixelList_m40(:,1)].*...
           double([RPixelValues_m40]))/(sum([RPixelValues_m40]));
   statsGwater(jj).WeightedCentroid_p40(2) = sum([PixelList_m40(:,2)].*...
           double([RPixelValues_m40]))/(sum([RPixelValues_m40]));
   statsGwater(jj).WeightedCentroid_p40(3) =  sum([PixelList_m40(:,3)].*...
           double([RPixelValues_m40]))/(sum([RPixelValues_m40]));
%also add these statxture, 
   centG_p40(jj,:) = statsGwater(jj).WeightedCentroid_p40;
   tintsG_p40(jj) = statsGwater(jj).tints_p40;
   volumeG_p40(jj) =statsGwater(jj).volume_p40;
   
   size2= 45;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m45 = [x,y,z];
   RPixelValues_m45 = zeros(numel(PixelList_m45(:,1)),1);
   for kk=1:numel(PixelList_m45(:,1))
      RPixelValues_m45(kk) = curr1b(PixelList_m45(kk,2),...
          PixelList_m45(kk,1), PixelList_m45(kk,3));
   end

   statsGwater(jj).tints_p45 = sum([RPixelValues_m45]);
   statsGwater(jj).area_p45 = numel([RPixelValues_m45]);
   statsGwater(jj).volume_p45 = numel([RPixelValues_m45]>0);
   statsGwater(jj).WeightedCentroid_p45(1) = sum([PixelList_m45(:,1)].*...
           double([RPixelValues_m45]))/(sum([RPixelValues_m45]));
   statsGwater(jj).WeightedCentroid_p45(2) = sum([PixelList_m45(:,2)].*...
           double([RPixelValues_m45]))/(sum([RPixelValues_m45]));
   statsGwater(jj).WeightedCentroid_p45(3) =  sum([PixelList_m45(:,3)].*...
           double([RPixelValues_m45]))/(sum([RPixelValues_m45]));
%also add these statxture, 
   centG_p45(jj,:) = statsGwater(jj).WeightedCentroid_p45;
   tintsG_p45(jj) = statsGwater(jj).tints_p45;
   volumeG_p45(jj) =statsGwater(jj).volume_p45;
   
   
   size2= 50;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m50 = [x,y,z];
   RPixelValues_m50 = zeros(numel(PixelList_m50(:,1)),1);
   for kk=1:numel(PixelList_m50(:,1))
      RPixelValues_m50(kk) = curr1b(PixelList_m50(kk,2),...
          PixelList_m50(kk,1), PixelList_m50(kk,3));
   end

   statsGwater(jj).tints_p50 = sum([RPixelValues_m50]);
   statsGwater(jj).area_p50 = numel([RPixelValues_m50]);
   statsGwater(jj).volume_p50 = numel([RPixelValues_m50]>0);
   statsGwater(jj).WeightedCentroid_p50(1) = sum([PixelList_m50(:,1)].*...
           double([RPixelValues_m50]))/(sum([RPixelValues_m50]));
   statsGwater(jj).WeightedCentroid_p50(2) = sum([PixelList_m50(:,2)].*...
           double([RPixelValues_m50]))/(sum([RPixelValues_m50]));
   statsGwater(jj).WeightedCentroid_p50(3) =  sum([PixelList_m50(:,3)].*...
           double([RPixelValues_m50]))/(sum([RPixelValues_m50]));
%also add these statxture, 
   centG_p50(jj,:) = statsGwater(jj).WeightedCentroid_p50;
   tintsG_p50(jj) = statsGwater(jj).tints_p50;
   volumeG_p50(jj) =statsGwater(jj).volume_p50;
   
   size2= 55;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m55 = [x,y,z];
   RPixelValues_m55 = zeros(numel(PixelList_m55(:,1)),1);
   for kk=1:numel(PixelList_m55(:,1))
      RPixelValues_m55(kk) = curr1b(PixelList_m55(kk,2),...
          PixelList_m55(kk,1), PixelList_m55(kk,3));
   end

   statsGwater(jj).tints_p55 = sum([RPixelValues_m55]);
   statsGwater(jj).area_p55 = numel([RPixelValues_m55]);
   statsGwater(jj).volume_p55 = numel([RPixelValues_m55]>0);
   statsGwater(jj).WeightedCentroid_p55(1) = sum([PixelList_m55(:,1)].*...
           double([RPixelValues_m55]))/(sum([RPixelValues_m55]));
   statsGwater(jj).WeightedCentroid_p55(2) = sum([PixelList_m55(:,2)].*...
           double([RPixelValues_m55]))/(sum([RPixelValues_m55]));
   statsGwater(jj).WeightedCentroid_p55(3) =  sum([PixelList_m55(:,3)].*...
           double([RPixelValues_m55]))/(sum([RPixelValues_m55]));
%also add these statxture, 
   centG_p55(jj,:) = statsGwater(jj).WeightedCentroid_p55;
   tintsG_p55(jj) = statsGwater(jj).tints_p55;
   volumeG_p55(jj) =statsGwater(jj).volume_p55;
   
   size2= 60;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m60 = [x,y,z];
   RPixelValues_m60 = zeros(numel(PixelList_m60(:,1)),1);
   for kk=1:numel(PixelList_m60(:,1))
      RPixelValues_m60(kk) = curr1b(PixelList_m60(kk,2),...
          PixelList_m60(kk,1), PixelList_m60(kk,3));
   end

   statsGwater(jj).tints_p60 = sum([RPixelValues_m60]);
   statsGwater(jj).area_p60 = numel([RPixelValues_m60]);
   statsGwater(jj).volume_p60 = numel([RPixelValues_m60]>0);
   statsGwater(jj).WeightedCentroid_p60(1) = sum([PixelList_m60(:,1)].*...
           double([RPixelValues_m60]))/(sum([RPixelValues_m60]));
   statsGwater(jj).WeightedCentroid_p60(2) = sum([PixelList_m60(:,2)].*...
           double([RPixelValues_m60]))/(sum([RPixelValues_m60]));
   statsGwater(jj).WeightedCentroid_p60(3) =  sum([PixelList_m60(:,3)].*...
           double([RPixelValues_m60]))/(sum([RPixelValues_m60]));
%also add these statxture, 
   centG_p60(jj,:) = statsGwater(jj).WeightedCentroid_p60;
   tintsG_p60(jj) = statsGwater(jj).tints_p60;
   volumeG_p60(jj) =statsGwater(jj).volume_p60;
   
   size2= 65;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m65 = [x,y,z];
   RPixelValues_m65 = zeros(numel(PixelList_m65(:,1)),1);
   for kk=1:numel(PixelList_m65(:,1))
      RPixelValues_m65(kk) = curr1b(PixelList_m65(kk,2),...
          PixelList_m65(kk,1), PixelList_m65(kk,3));
   end

   statsGwater(jj).tints_p65 = sum([RPixelValues_m65]);
   statsGwater(jj).area_p65 = numel([RPixelValues_m65]);
   statsGwater(jj).volume_p65 = numel([RPixelValues_m65]>0);
   statsGwater(jj).WeightedCentroid_p65(1) = sum([PixelList_m65(:,1)].*...
           double([RPixelValues_m65]))/(sum([RPixelValues_m65]));
   statsGwater(jj).WeightedCentroid_p65(2) = sum([PixelList_m65(:,2)].*...
           double([RPixelValues_m65]))/(sum([RPixelValues_m65]));
   statsGwater(jj).WeightedCentroid_p65(3) =  sum([PixelList_m65(:,3)].*...
           double([RPixelValues_m65]))/(sum([RPixelValues_m65]));
%also add these statxture, 
   centG_p65(jj,:) = statsGwater(jj).WeightedCentroid_p65;
   tintsG_p65(jj) = statsGwater(jj).tints_p65;
   volumeG_p65(jj) =statsGwater(jj).volume_p65;

   size2= 70;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m70 = [x,y,z];
   RPixelValues_m70 = zeros(numel(PixelList_m70(:,1)),1);
   for kk=1:numel(PixelList_m70(:,1))
      RPixelValues_m70(kk) = curr1b(PixelList_m70(kk,2),...
          PixelList_m70(kk,1), PixelList_m70(kk,3));
   end

   statsGwater(jj).tints_p70 = sum([RPixelValues_m70]);
   statsGwater(jj).area_p70 = numel([RPixelValues_m70]);
   statsGwater(jj).volume_p70 = numel([RPixelValues_m70]>0);
   statsGwater(jj).WeightedCentroid_p70(1) = sum([PixelList_m70(:,1)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   statsGwater(jj).WeightedCentroid_p70(2) = sum([PixelList_m70(:,2)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   statsGwater(jj).WeightedCentroid_p70(3) =  sum([PixelList_m70(:,3)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
%also add these statxture, 
   centG_p70(jj,:) = statsGwater(jj).WeightedCentroid_p70;
   tintsG_p70(jj) = statsGwater(jj).tints_p70;
   volumeG_p70(jj) =statsGwater(jj).volume_p70;

end
%
Result = zeros(numel(statsGwater),1);
for i = 1:numel(statsGwater)
    x = 1:1:14;
    x = x';
    a = zeros(14,1);
    a(1) = statsGwater(i).area_p5 / statsGwater(i).area_p70;
    a(2) = statsGwater(i).area_p10 / statsGwater(i).area_p70;
    a(3) = statsGwater(i).area_p15 / statsGwater(i).area_p70;
    a(4) = statsGwater(i).area_p20 / statsGwater(i).area_p70;
    a(5) = statsGwater(i).area_p25 / statsGwater(i).area_p70;
    a(6) = statsGwater(i).area_p30 / statsGwater(i).area_p70;
    a(7) = statsGwater(i).area_p35 / statsGwater(i).area_p70;
    a(8) = statsGwater(i).area_p40 / statsGwater(i).area_p70;
    a(9) = statsGwater(i).area_p45 / statsGwater(i).area_p70;
    a(10) = statsGwater(i).area_p50 / statsGwater(i).area_p70;
    a(11) = statsGwater(i).area_p55 / statsGwater(i).area_p70;
    a(12) = statsGwater(i).area_p60 / statsGwater(i).area_p70;
    a(13) = statsGwater(i).area_p65 / statsGwater(i).area_p70;
    a(14) = statsGwater(i).area_p70 / statsGwater(i).area_p70;
    [f,gof] = fit(x,a,'smoothingspline');
    R = feval(f,1:0.01:14);
    x = find(R>0.5);
    Result(i) = 1 + x(1) * 0.01;
end
Result_G = Result;

%
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

outpath = [exp_folder 'Result\'];
outpath_Vsyn = [outpath '5_V_Syn\'];
%
load([outpath_Vsyn 'G_paired_3.mat']);
statsRwater = statsGwater_ssss;
clear statsGwater
%
clear statsGa2s statsGwater
load([outpath_Vsyn 'R_paired_3.mat']);
statsGwater = statsRwater_ssss;
%
% num_images = num_images - 2;
clear BP BG bg2
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
BP2 = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsGwater)
    disp(int2str(i))
    BP(statsGwater(i).PixelIdxList)=statsGwater(i).PixelValues;  
end
%
%BP2 = BP;
for i = 1:numel(statsRwater)
    disp(int2str(i))
    BP2(statsRwater(i).PixelIdxList)=statsRwater(i).PixelValues;  
end
%
disp('loading data')
%
%if not present, load in stats list

%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsGwater):-1:1
    disp(i)
    
statsGwater(i).tints_p5 = [];
statsGwater(i).volume_p5 = [];
statsGwater(i).area_p5 = [];
statsGwater(i).WeightedCentroid_p5 = [];
centG_p5(i,1:3) = 0;
tintsG_p5(i) = 0;
volumeG_p5(i) =0;

statsGwater(i).tints_p10 = [];
statsGwater(i).volume_p10 = [];
statsGwater(i).area_p10 = [];
statsGwater(i).WeightedCentroid_p10 = [];
centG_p10(i,1:3) = 0;
tintsG_p10(i) = 0;
volumeG_p10(i) =0;

statsGwater(i).tints_p15 = [];
statsGwater(i).volume_p15 = [];
statsGwater(i).area_p15 = [];
statsGwater(i).WeightedCentroid_p15 = [];
centG_p15(i,1:3) = 0;
tintsG_p15(i) = 0;
volumeG_p15(i) =0;

statsGwater(i).tints_p20 = [];
statsGwater(i).volume_p20 = [];
statsGwater(i).area_p20 = [];
statsGwater(i).WeightedCentroid_p20 = [];
centG_p20(i,1:3) = 0;
tintsG_p20(i) = 0;
volumeG_p20(i) =0;

statsGwater(i).tints_p25 = [];
statsGwater(i).volume_p25 = [];
statsGwater(i).area_p25 = [];
statsGwater(i).WeightedCentroid_p25 = [];
centG_p25(i,1:3) = 0;
tintsG_p25(i) = 0;
volumeG_p25(i) =0;


statsGwater(i).tints_p30 = [];
statsGwater(i).volume_p30 = [];
statsGwater(i).area_p30 = [];
statsGwater(i).WeightedCentroid_p30 = [];
centG_p30(i,1:3) = 0;
tintsG_p30(i) = 0;
volumeG_p30(i) =0;

statsGwater(i).tints_p35 = [];
statsGwater(i).volume_p35 = [];
statsGwater(i).area_p35 = [];
statsGwater(i).WeightedCentroid_p35 = [];
centG_p35(i,1:3) = 0;
tintsG_p35(i) = 0;
volumeG_p35(i) =0;

statsGwater(i).tints_p40 = [];
statsGwater(i).volume_p40 = [];
statsGwater(i).area_p40 = [];
statsGwater(i).WeightedCentroid_p40 = [];
centG_p40(i,1:3) = 0;
tintsG_p40(i) = 0;
volumeG_p40(i) =0;

statsGwater(i).tints_p45 = [];
statsGwater(i).volume_p45 = [];
statsGwater(i).area_p45 = [];
statsGwater(i).WeightedCentroid_p45 = [];
centG_p45(i,1:3) = 0;
tintsG_p45(i) = 0;
volumeG_p45(i) =0;

statsGwater(i).tints_p50 = [];
statsGwater(i).volume_p50 = [];
statsGwater(i).area_p50 = [];
statsGwater(i).WeightedCentroid_p50 = [];
centG_p50(i,1:3) = 0;
tintsG_p50(i) = 0;
volumeG_p50(i) =0;

statsGwater(i).tints_p55 = [];
statsGwater(i).volume_p55 = [];
statsGwater(i).area_p55 = [];
statsGwater(i).WeightedCentroid_p55 = [];
centG_p55(i,1:3) = 0;
tintsG_p55(i) = 0;
volumeG_p55(i) =0;

statsGwater(i).tints_p60 = [];
statsGwater(i).volume_p60 = [];
statsGwater(i).area_p60 = [];
statsGwater(i).WeightedCentroid_p60 = [];
centG_p60(i,1:3) = 0;
tintsG_p60(i) = 0;
volumeG_p60(i) =0;

statsGwater(i).tints_p65 = [];
statsGwater(i).volume_p65 = [];
statsGwater(i).area_p65 = [];
statsGwater(i).WeightedCentroid_p65 = [];
centG_p65(i,1:3) = 0;
tintsG_p65(i) = 0;
volumeG_p65(i) =0;

statsGwater(i).tints_p70 = [];
statsGwater(i).volume_p70 = [];
statsGwater(i).area_p70 = [];
statsGwater(i).WeightedCentroid_p70 = [];
centG_p70(i,1:3) = 0;
tintsG_p70(i) = 0;
volumeG_p70(i) =0;
%try preproccessing the BP varialbe in a for loop
   minpix = min(statsGwater(i).PixelList);  maxpix = max(statsGwater(i).PixelList);
   min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
   max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end      
   BPp(i).mat = BP(min2:max2,min1:max1,min3:max3); 
   BP2p(i).mat = BP2(min2:max2,min1:max1,min3:max3); 
end
%

for jj=1:numel(statsGwater)
   %
   disp(jj)
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
   curr1a = BPp(jj).mat;
   curr1b = BP2p(jj).mat;
   Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
%
    size2= 5;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m5 = [x,y,z];
   RPixelValues_m5 = zeros(numel(PixelList_m5(:,1)),1);
   for kk=1:numel(PixelList_m5(:,1))
      RPixelValues_m5(kk) = curr1b(PixelList_m5(kk,2),...
          PixelList_m5(kk,1), PixelList_m5(kk,3));
   end

   statsGwater(jj).tints_p5 = sum([RPixelValues_m5]);
   statsGwater(jj).area_p5 = numel([RPixelValues_m5]);
   statsGwater(jj).volume_p5 = numel([RPixelValues_m5]>0);
   statsGwater(jj).WeightedCentroid_p5(1) = sum([PixelList_m5(:,1)].*...
           double([RPixelValues_m5]))/(sum([RPixelValues_m5]));
   statsGwater(jj).WeightedCentroid_p5(2) = sum([PixelList_m5(:,2)].*...
           double([RPixelValues_m5]))/(sum([RPixelValues_m5]));
   statsGwater(jj).WeightedCentroid_p5(3) =  sum([PixelList_m5(:,3)].*...
           double([RPixelValues_m5]))/(sum([RPixelValues_m5]));
%also add these statxture, 
   centG_p5(jj,:) = statsGwater(jj).WeightedCentroid_p5;
   tintsG_p5(jj) = statsGwater(jj).tints_p5;
   volumeG_p5(jj) =statsGwater(jj).volume_p5;
   
   
   size2= 10;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m10 = [x,y,z];
   RPixelValues_m10 = zeros(numel(PixelList_m10(:,1)),1);
   for kk=1:numel(PixelList_m10(:,1))
      RPixelValues_m10(kk) = curr1b(PixelList_m10(kk,2),...
          PixelList_m10(kk,1), PixelList_m10(kk,3));
   end

   statsGwater(jj).tints_p10 = sum([RPixelValues_m10]);
   statsGwater(jj).area_p10 = numel([RPixelValues_m10]);
   statsGwater(jj).volume_p10 = numel([RPixelValues_m10]>0);
   statsGwater(jj).WeightedCentroid_p10(1) = sum([PixelList_m10(:,1)].*...
           double([RPixelValues_m10]))/(sum([RPixelValues_m10]));
   statsGwater(jj).WeightedCentroid_p10(2) = sum([PixelList_m10(:,2)].*...
           double([RPixelValues_m10]))/(sum([RPixelValues_m10]));
   statsGwater(jj).WeightedCentroid_p10(3) =  sum([PixelList_m10(:,3)].*...
           double([RPixelValues_m10]))/(sum([RPixelValues_m10]));
%also add these statxture, 
   centG_p10(jj,:) = statsGwater(jj).WeightedCentroid_p10;
   tintsG_p10(jj) = statsGwater(jj).tints_p10;
   volumeG_p10(jj) =statsGwater(jj).volume_p10;

   size2= 15;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m15 = [x,y,z];
   RPixelValues_m15 = zeros(numel(PixelList_m15(:,1)),1);
   for kk=1:numel(PixelList_m15(:,1))
      RPixelValues_m15(kk) = curr1b(PixelList_m15(kk,2),...
          PixelList_m15(kk,1), PixelList_m15(kk,3));
   end

   statsGwater(jj).tints_p15 = sum([RPixelValues_m15]);
   statsGwater(jj).area_p15 = numel([RPixelValues_m15]);
   statsGwater(jj).volume_p15 = numel([RPixelValues_m15]>0);
   statsGwater(jj).WeightedCentroid_p15(1) = sum([PixelList_m15(:,1)].*...
           double([RPixelValues_m15]))/(sum([RPixelValues_m15]));
   statsGwater(jj).WeightedCentroid_p15(2) = sum([PixelList_m15(:,2)].*...
           double([RPixelValues_m15]))/(sum([RPixelValues_m15]));
   statsGwater(jj).WeightedCentroid_p15(3) =  sum([PixelList_m15(:,3)].*...
           double([RPixelValues_m15]))/(sum([RPixelValues_m15]));
%also add these statxture, 
   centG_p15(jj,:) = statsGwater(jj).WeightedCentroid_p15;
   tintsG_p15(jj) = statsGwater(jj).tints_p15;
   volumeG_p15(jj) =statsGwater(jj).volume_p15;
   
   size2= 20;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m20 = [x,y,z];
   RPixelValues_m20 = zeros(numel(PixelList_m20(:,1)),1);
   for kk=1:numel(PixelList_m20(:,1))
      RPixelValues_m20(kk) = curr1b(PixelList_m20(kk,2),...
          PixelList_m20(kk,1), PixelList_m20(kk,3));
   end

   statsGwater(jj).tints_p20 = sum([RPixelValues_m20]);
   statsGwater(jj).area_p20 = numel([RPixelValues_m20]);
   statsGwater(jj).volume_p20 = numel([RPixelValues_m20]>0);
   statsGwater(jj).WeightedCentroid_p20(1) = sum([PixelList_m20(:,1)].*...
           double([RPixelValues_m20]))/(sum([RPixelValues_m20]));
   statsGwater(jj).WeightedCentroid_p20(2) = sum([PixelList_m20(:,2)].*...
           double([RPixelValues_m20]))/(sum([RPixelValues_m20]));
   statsGwater(jj).WeightedCentroid_p20(3) =  sum([PixelList_m20(:,3)].*...
           double([RPixelValues_m20]))/(sum([RPixelValues_m20]));
%also add these statxture, 
   centG_p20(jj,:) = statsGwater(jj).WeightedCentroid_p20;
   tintsG_p20(jj) = statsGwater(jj).tints_p20;
   volumeG_p20(jj) =statsGwater(jj).volume_p20;
   
   size2= 25;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m25 = [x,y,z];
   RPixelValues_m25 = zeros(numel(PixelList_m25(:,1)),1);
   for kk=1:numel(PixelList_m25(:,1))
      RPixelValues_m25(kk) = curr1b(PixelList_m25(kk,2),...
          PixelList_m25(kk,1), PixelList_m25(kk,3));
   end

   statsGwater(jj).tints_p25 = sum([RPixelValues_m25]);
   statsGwater(jj).area_p25 = numel([RPixelValues_m25]);
   statsGwater(jj).volume_p25 = numel([RPixelValues_m25]>0);
   statsGwater(jj).WeightedCentroid_p25(1) = sum([PixelList_m25(:,1)].*...
           double([RPixelValues_m25]))/(sum([RPixelValues_m25]));
   statsGwater(jj).WeightedCentroid_p25(2) = sum([PixelList_m25(:,2)].*...
           double([RPixelValues_m25]))/(sum([RPixelValues_m25]));
   statsGwater(jj).WeightedCentroid_p25(3) =  sum([PixelList_m25(:,3)].*...
           double([RPixelValues_m25]))/(sum([RPixelValues_m25]));
%also add these statxture, 
   centG_p25(jj,:) = statsGwater(jj).WeightedCentroid_p25;
   tintsG_p25(jj) = statsGwater(jj).tints_p25;
   volumeG_p25(jj) =statsGwater(jj).volume_p25;
   
   
   size2= 30;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m30 = [x,y,z];
   RPixelValues_m30 = zeros(numel(PixelList_m30(:,1)),1);
   for kk=1:numel(PixelList_m30(:,1))
      RPixelValues_m30(kk) = curr1b(PixelList_m30(kk,2),...
          PixelList_m30(kk,1), PixelList_m30(kk,3));
   end

   statsGwater(jj).tints_p30 = sum([RPixelValues_m30]);
   statsGwater(jj).area_p30 = numel([RPixelValues_m30]);
   statsGwater(jj).volume_p30 = numel([RPixelValues_m30]>0);
   statsGwater(jj).WeightedCentroid_p30(1) = sum([PixelList_m30(:,1)].*...
           double([RPixelValues_m30]))/(sum([RPixelValues_m30]));
   statsGwater(jj).WeightedCentroid_p30(2) = sum([PixelList_m30(:,2)].*...
           double([RPixelValues_m30]))/(sum([RPixelValues_m30]));
   statsGwater(jj).WeightedCentroid_p30(3) =  sum([PixelList_m30(:,3)].*...
           double([RPixelValues_m30]))/(sum([RPixelValues_m30]));
%also add these statxture, 
   centG_p30(jj,:) = statsGwater(jj).WeightedCentroid_p30;
   tintsG_p30(jj) = statsGwater(jj).tints_p30;
   volumeG_p30(jj) =statsGwater(jj).volume_p30;
   
   size2= 35;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m35 = [x,y,z];
   RPixelValues_m35 = zeros(numel(PixelList_m35(:,1)),1);
   for kk=1:numel(PixelList_m35(:,1))
      RPixelValues_m35(kk) = curr1b(PixelList_m35(kk,2),...
          PixelList_m35(kk,1), PixelList_m35(kk,3));
   end

   statsGwater(jj).tints_p35 = sum([RPixelValues_m35]);
   statsGwater(jj).area_p35 = numel([RPixelValues_m35]);
   statsGwater(jj).volume_p35 = numel([RPixelValues_m35]>0);
   statsGwater(jj).WeightedCentroid_p35(1) = sum([PixelList_m35(:,1)].*...
           double([RPixelValues_m35]))/(sum([RPixelValues_m35]));
   statsGwater(jj).WeightedCentroid_p35(2) = sum([PixelList_m35(:,2)].*...
           double([RPixelValues_m35]))/(sum([RPixelValues_m35]));
   statsGwater(jj).WeightedCentroid_p35(3) =  sum([PixelList_m35(:,3)].*...
           double([RPixelValues_m35]))/(sum([RPixelValues_m35]));
%also add these statxture, 
   centG_p35(jj,:) = statsGwater(jj).WeightedCentroid_p35;
   tintsG_p35(jj) = statsGwater(jj).tints_p35;
   volumeG_p35(jj) =statsGwater(jj).volume_p35;
   
   
   size2= 40;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m40 = [x,y,z];
   RPixelValues_m40 = zeros(numel(PixelList_m40(:,1)),1);
   for kk=1:numel(PixelList_m40(:,1))
      RPixelValues_m40(kk) = curr1b(PixelList_m40(kk,2),...
          PixelList_m40(kk,1), PixelList_m40(kk,3));
   end

   statsGwater(jj).tints_p40 = sum([RPixelValues_m40]);
   statsGwater(jj).area_p40 = numel([RPixelValues_m40]);
   statsGwater(jj).volume_p40 = numel([RPixelValues_m40]>0);
   statsGwater(jj).WeightedCentroid_p40(1) = sum([PixelList_m40(:,1)].*...
           double([RPixelValues_m40]))/(sum([RPixelValues_m40]));
   statsGwater(jj).WeightedCentroid_p40(2) = sum([PixelList_m40(:,2)].*...
           double([RPixelValues_m40]))/(sum([RPixelValues_m40]));
   statsGwater(jj).WeightedCentroid_p40(3) =  sum([PixelList_m40(:,3)].*...
           double([RPixelValues_m40]))/(sum([RPixelValues_m40]));
%also add these statxture, 
   centG_p40(jj,:) = statsGwater(jj).WeightedCentroid_p40;
   tintsG_p40(jj) = statsGwater(jj).tints_p40;
   volumeG_p40(jj) =statsGwater(jj).volume_p40;
   
   size2= 45;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m45 = [x,y,z];
   RPixelValues_m45 = zeros(numel(PixelList_m45(:,1)),1);
   for kk=1:numel(PixelList_m45(:,1))
      RPixelValues_m45(kk) = curr1b(PixelList_m45(kk,2),...
          PixelList_m45(kk,1), PixelList_m45(kk,3));
   end

   statsGwater(jj).tints_p45 = sum([RPixelValues_m45]);
   statsGwater(jj).area_p45 = numel([RPixelValues_m45]);
   statsGwater(jj).volume_p45 = numel([RPixelValues_m45]>0);
   statsGwater(jj).WeightedCentroid_p45(1) = sum([PixelList_m45(:,1)].*...
           double([RPixelValues_m45]))/(sum([RPixelValues_m45]));
   statsGwater(jj).WeightedCentroid_p45(2) = sum([PixelList_m45(:,2)].*...
           double([RPixelValues_m45]))/(sum([RPixelValues_m45]));
   statsGwater(jj).WeightedCentroid_p45(3) =  sum([PixelList_m45(:,3)].*...
           double([RPixelValues_m45]))/(sum([RPixelValues_m45]));
%also add these statxture, 
   centG_p45(jj,:) = statsGwater(jj).WeightedCentroid_p45;
   tintsG_p45(jj) = statsGwater(jj).tints_p45;
   volumeG_p45(jj) =statsGwater(jj).volume_p45;
   
   
   size2= 50;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m50 = [x,y,z];
   RPixelValues_m50 = zeros(numel(PixelList_m50(:,1)),1);
   for kk=1:numel(PixelList_m50(:,1))
      RPixelValues_m50(kk) = curr1b(PixelList_m50(kk,2),...
          PixelList_m50(kk,1), PixelList_m50(kk,3));
   end

   statsGwater(jj).tints_p50 = sum([RPixelValues_m50]);
   statsGwater(jj).area_p50 = numel([RPixelValues_m50]);
   statsGwater(jj).volume_p50 = numel([RPixelValues_m50]>0);
   statsGwater(jj).WeightedCentroid_p50(1) = sum([PixelList_m50(:,1)].*...
           double([RPixelValues_m50]))/(sum([RPixelValues_m50]));
   statsGwater(jj).WeightedCentroid_p50(2) = sum([PixelList_m50(:,2)].*...
           double([RPixelValues_m50]))/(sum([RPixelValues_m50]));
   statsGwater(jj).WeightedCentroid_p50(3) =  sum([PixelList_m50(:,3)].*...
           double([RPixelValues_m50]))/(sum([RPixelValues_m50]));
%also add these statxture, 
   centG_p50(jj,:) = statsGwater(jj).WeightedCentroid_p50;
   tintsG_p50(jj) = statsGwater(jj).tints_p50;
   volumeG_p50(jj) =statsGwater(jj).volume_p50;
   
   size2= 55;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m55 = [x,y,z];
   RPixelValues_m55 = zeros(numel(PixelList_m55(:,1)),1);
   for kk=1:numel(PixelList_m55(:,1))
      RPixelValues_m55(kk) = curr1b(PixelList_m55(kk,2),...
          PixelList_m55(kk,1), PixelList_m55(kk,3));
   end

   statsGwater(jj).tints_p55 = sum([RPixelValues_m55]);
   statsGwater(jj).area_p55 = numel([RPixelValues_m55]);
   statsGwater(jj).volume_p55 = numel([RPixelValues_m55]>0);
   statsGwater(jj).WeightedCentroid_p55(1) = sum([PixelList_m55(:,1)].*...
           double([RPixelValues_m55]))/(sum([RPixelValues_m55]));
   statsGwater(jj).WeightedCentroid_p55(2) = sum([PixelList_m55(:,2)].*...
           double([RPixelValues_m55]))/(sum([RPixelValues_m55]));
   statsGwater(jj).WeightedCentroid_p55(3) =  sum([PixelList_m55(:,3)].*...
           double([RPixelValues_m55]))/(sum([RPixelValues_m55]));
%also add these statxture, 
   centG_p55(jj,:) = statsGwater(jj).WeightedCentroid_p55;
   tintsG_p55(jj) = statsGwater(jj).tints_p55;
   volumeG_p55(jj) =statsGwater(jj).volume_p55;
   
   size2= 60;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m60 = [x,y,z];
   RPixelValues_m60 = zeros(numel(PixelList_m60(:,1)),1);
   for kk=1:numel(PixelList_m60(:,1))
      RPixelValues_m60(kk) = curr1b(PixelList_m60(kk,2),...
          PixelList_m60(kk,1), PixelList_m60(kk,3));
   end

   statsGwater(jj).tints_p60 = sum([RPixelValues_m60]);
   statsGwater(jj).area_p60 = numel([RPixelValues_m60]);
   statsGwater(jj).volume_p60 = numel([RPixelValues_m60]>0);
   statsGwater(jj).WeightedCentroid_p60(1) = sum([PixelList_m60(:,1)].*...
           double([RPixelValues_m60]))/(sum([RPixelValues_m60]));
   statsGwater(jj).WeightedCentroid_p60(2) = sum([PixelList_m60(:,2)].*...
           double([RPixelValues_m60]))/(sum([RPixelValues_m60]));
   statsGwater(jj).WeightedCentroid_p60(3) =  sum([PixelList_m60(:,3)].*...
           double([RPixelValues_m60]))/(sum([RPixelValues_m60]));
%also add these statxture, 
   centG_p60(jj,:) = statsGwater(jj).WeightedCentroid_p60;
   tintsG_p60(jj) = statsGwater(jj).tints_p60;
   volumeG_p60(jj) =statsGwater(jj).volume_p60;
   
   size2= 65;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m65 = [x,y,z];
   RPixelValues_m65 = zeros(numel(PixelList_m65(:,1)),1);
   for kk=1:numel(PixelList_m65(:,1))
      RPixelValues_m65(kk) = curr1b(PixelList_m65(kk,2),...
          PixelList_m65(kk,1), PixelList_m65(kk,3));
   end

   statsGwater(jj).tints_p65 = sum([RPixelValues_m65]);
   statsGwater(jj).area_p65 = numel([RPixelValues_m65]);
   statsGwater(jj).volume_p65 = numel([RPixelValues_m65]>0);
   statsGwater(jj).WeightedCentroid_p65(1) = sum([PixelList_m65(:,1)].*...
           double([RPixelValues_m65]))/(sum([RPixelValues_m65]));
   statsGwater(jj).WeightedCentroid_p65(2) = sum([PixelList_m65(:,2)].*...
           double([RPixelValues_m65]))/(sum([RPixelValues_m65]));
   statsGwater(jj).WeightedCentroid_p65(3) =  sum([PixelList_m65(:,3)].*...
           double([RPixelValues_m65]))/(sum([RPixelValues_m65]));
%also add these statxture, 
   centG_p65(jj,:) = statsGwater(jj).WeightedCentroid_p65;
   tintsG_p65(jj) = statsGwater(jj).tints_p65;
   volumeG_p65(jj) =statsGwater(jj).volume_p65;

   size2= 70;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m70 = [x,y,z];
   RPixelValues_m70 = zeros(numel(PixelList_m70(:,1)),1);
   for kk=1:numel(PixelList_m70(:,1))
      RPixelValues_m70(kk) = curr1b(PixelList_m70(kk,2),...
          PixelList_m70(kk,1), PixelList_m70(kk,3));
   end

   statsGwater(jj).tints_p70 = sum([RPixelValues_m70]);
   statsGwater(jj).area_p70 = numel([RPixelValues_m70]);
   statsGwater(jj).volume_p70 = numel([RPixelValues_m70]>0);
   statsGwater(jj).WeightedCentroid_p70(1) = sum([PixelList_m70(:,1)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   statsGwater(jj).WeightedCentroid_p70(2) = sum([PixelList_m70(:,2)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   statsGwater(jj).WeightedCentroid_p70(3) =  sum([PixelList_m70(:,3)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
%also add these statxture, 
   centG_p70(jj,:) = statsGwater(jj).WeightedCentroid_p70;
   tintsG_p70(jj) = statsGwater(jj).tints_p70;
   volumeG_p70(jj) =statsGwater(jj).volume_p70;

end
%
Result = zeros(numel(statsGwater),1);
for i = 1:numel(statsGwater)
    x = 1:1:14;
    x = x';
    a = zeros(14,1);
    a(1) = statsGwater(i).area_p5 / statsGwater(i).area_p70;
    a(2) = statsGwater(i).area_p10 / statsGwater(i).area_p70;
    a(3) = statsGwater(i).area_p15 / statsGwater(i).area_p70;
    a(4) = statsGwater(i).area_p20 / statsGwater(i).area_p70;
    a(5) = statsGwater(i).area_p25 / statsGwater(i).area_p70;
    a(6) = statsGwater(i).area_p30 / statsGwater(i).area_p70;
    a(7) = statsGwater(i).area_p35 / statsGwater(i).area_p70;
    a(8) = statsGwater(i).area_p40 / statsGwater(i).area_p70;
    a(9) = statsGwater(i).area_p45 / statsGwater(i).area_p70;
    a(10) = statsGwater(i).area_p50 / statsGwater(i).area_p70;
    a(11) = statsGwater(i).area_p55 / statsGwater(i).area_p70;
    a(12) = statsGwater(i).area_p60 / statsGwater(i).area_p70;
    a(13) = statsGwater(i).area_p65 / statsGwater(i).area_p70;
    a(14) = statsGwater(i).area_p70 / statsGwater(i).area_p70;
    [f,gof] = fit(x,a,'smoothingspline');
    R = feval(f,1:0.01:14);
    x = find(R>0.5);
    Result(i) = 1 + x(1) * 0.01;
end
Result_R = Result;
%
save([outpath '6_Dist/Add_To_Dist.mat'],'Result_G','Result_R');