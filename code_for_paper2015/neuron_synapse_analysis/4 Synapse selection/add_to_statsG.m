%add_to_statsG
%add statxture and pixellist_140 and opposing channel volume, tints, and 
%weighted centroid to statslist
%load in statslist and opposing channel image
%cluster_seg_whole
%%
clear all
%%
base_path = '/n/contefs1/backup/ysigal/';
exp_folder = [base_path 'GLD3_bistratified_RGC/'];
%path = [exp_folder  'proc_2nd_paper/cropped/'];
%mergedpath = [path 'ipl_synaptic_pre2/sections/'];
%mergedpath2 = [path 'storm_merged_thresh/'];
path = [exp_folder  'elastic_thresh/'];
mergedpath = [path 'storm_merged/'];
%
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);
disp(num_images)
%set voxel size in nm
voxel=[15.8, 15.8, 70];
%%
%
load([path 'statsP2sw10_edges.mat']);
load([path 'statsP2nsw10_edges.mat']);
%%
load([path 'statsG2sw10.mat']);
%%
parpool(32)
%%
clear BP BG bg2
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsPa1s)
    disp(int2str(i))
    BP(statsPa1s(i).PixelIdxList)=statsPa1s(i).PixelValues;  
end
BP2 = BP;
for i = 1:numel(statsPa1ns)
    disp(int2str(i))
    BP2(statsPa1ns(i).PixelIdxList)=statsPa1ns(i).PixelValues;  
end
disp('loading data')

%if not present, load in stats list

statsGwater = statsGa2s;
%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsGwater):-1:1
    disp(i)
statsGwater(i).tints_p20 = [];
statsGwater(i).volume_p20 = [];
statsGwater(i).area_p20 = [];
statsGwater(i).WeightedCentroid_p20 = [];
statsGwater(i).statxtureP_20all = [];
statsGwater(i).statxtureP_20pos = [];
centG_p20(i,1:3) = 0;
tintsG_p20(i) = 0;
volumeG_p20(i) =0;

statsGwater(i).tints_p60 = [];
statsGwater(i).volume_p60 = [];
statsGwater(i).area_p60 = [];
statsGwater(i).WeightedCentroid_p60 = [];
statsGwater(i).statxtureP_60all = [];
statsGwater(i).statxtureP_60pos = [];
centG_p60(i,1:3) = 0;
tintsG_p60(i) = 0;
volumeG_p60(i) =0;

statsGwater(i).tints_p100 = [];
statsGwater(i).volume_p100 = [];
statsGwater(i).area_p100 = [];
statsGwater(i).WeightedCentroid_p100 = [];
statsGwater(i).statxtureP_100all = [];
statsGwater(i).statxtureP_100pos = [];
centG_p100(i,1:3) = 0;
tintsG_p100(i) = 0;
volumeG_p100(i) =0;

statsGwater(i).tints_p140 = [];
statsGwater(i).volume_p140 = [];
statsGwater(i).area_p140 = [];
statsGwater(i).WeightedCentroid_p140 = [];
statsGwater(i).statxtureP_140all = [];
statsGwater(i).statxtureP_140pos = [];
centG_p140(i,1:3) = 0;
tintsG_p140(i) = 0;
volumeG_p140(i) =0;

statsGwater(i).tints_P140 = [];
statsGwater(i).volume_P140 = [];
statsGwater(i).area_P140 = [];
statsGwater(i).WeightedCentroid_P140 = [];
statsGwater(i).statxtureP_P140all = [];
statsGwater(i).statxtureP_P140pos = [];
centG_P140(i,1:3) = 0;
tintsG_P140(i) = 0;
volumeG_P140(i) =0;

statsGwater(i).tints_p180 = [];
statsGwater(i).volume_p180 = [];
statsGwater(i).area_p180 = [];
statsGwater(i).WeightedCentroid_p180 = [];
statsGwater(i).statxtureP_180all = [];
statsGwater(i).statxtureP_180pos = [];
centG_p180(i,1:3) = 0;
tintsG_p180(i) = 0;
volumeG_p180(i) =0;

statsGwater(i).tints_p220 = [];
statsGwater(i).volume_p220 = [];
statsGwater(i).area_p220 = [];
statsGwater(i).WeightedCentroid_p220 = [];
statsGwater(i).statxtureP_220all = [];
statsGwater(i).statxtureP_220pos = [];
centG_p220(i,1:3) = 0;
tintsG_p220(i) = 0;
volumeG_p220(i) =0;

statsGwater(i).tints_p260 = [];
statsGwater(i).volume_p260 = [];
statsGwater(i).area_p260 = [];
statsGwater(i).WeightedCentroid_p260 = [];
statsGwater(i).statxtureP_260all = [];
statsGwater(i).statxtureP_260pos = [];
centG_p260(i,1:3) = 0;
tintsG_p260(i) = 0;
volumeG_p260(i) =0;

statsGwater(i).tints_p300 = [];
statsGwater(i).volume_p300 = [];
statsGwater(i).area_p300 = [];
statsGwater(i).WeightedCentroid_p300 = [];
statsGwater(i).statxtureP_300all = [];
statsGwater(i).statxtureP_300pos = [];
centG_p300(i,1:3) = 0;
tintsG_p300(i) = 0;
volumeG_p300(i) =0;

statsGwater(i).tints_p400 = [];
statsGwater(i).volume_p400 = [];
statsGwater(i).area_p400 = [];
statsGwater(i).WeightedCentroid_p400 = [];
statsGwater(i).statxtureP_400all = [];
statsGwater(i).statxtureP_400pos = [];
centG_p400(i,1:3) = 0;
tintsG_p400(i) = 0;
volumeG_p400(i) =0;

statsGwater(i).tints_p500 = [];
statsGwater(i).volume_p500 = [];
statsGwater(i).area_p500 = [];
statsGwater(i).WeightedCentroid_p500 = [];
statsGwater(i).statxtureP_500all = [];
statsGwater(i).statxtureP_500pos = [];
centG_p500(i,1:3) = 0;
tintsG_p500(i) = 0;
volumeG_p500(i) =0;


statsGwater(i).statxture_all = [];
statsGwater(i).statxture_pos = [];
statxture_all(i,1:6) = 0;
statxture_pos(i,1:6) = 0;

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
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
parfor jj=1:numel(statsGwater)
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

   size2= 20;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m20 = [x,y,z];
   RPixelValues_m20 = zeros(numel(PixelList_m20(:,1)),1);
   for kk=1:numel(PixelList_m20(:,1))
      RPixelValues_m20(kk) = curr1a(PixelList_m20(kk,2),...
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
   statsGwater(jj).statxtureP_20all = statxture2([RPixelValues_m20]);
   posvals2 = [RPixelValues_m20]>0;
   statsGwater(jj).statxtureP_20pos = statxture2([RPixelValues_m20(posvals2)]);
   centG_p20(jj,:) = statsGwater(jj).WeightedCentroid_p20;
   tintsG_p20(jj) = statsGwater(jj).tints_p20;
   volumeG_p20(jj) =statsGwater(jj).volume_p20;
    
   size2= 60;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m60 = [x,y,z];
   RPixelValues_m60 = zeros(numel(PixelList_m60(:,1)),1);
   for kk=1:numel(PixelList_m60(:,1))
      RPixelValues_m60(kk) = curr1a(PixelList_m60(kk,2),...
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
   statsGwater(jj).statxtureP_60all = statxture2([RPixelValues_m60]);
   posvals2 = [RPixelValues_m60]>0;
   statsGwater(jj).statxtureP_60pos = statxture2([RPixelValues_m60(posvals2)]);
   centG_p60(jj,:) = statsGwater(jj).WeightedCentroid_p60;
   tintsG_p60(jj) = statsGwater(jj).tints_p60;
   volumeG_p60(jj) =statsGwater(jj).volume_p60;
    
   size2= 100;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m100 = [x,y,z];
   RPixelValues_m100 = zeros(numel(PixelList_m100(:,1)),1);
   for kk=1:numel(PixelList_m100(:,1))
      RPixelValues_m100(kk) = curr1a(PixelList_m100(kk,2),...
          PixelList_m100(kk,1), PixelList_m100(kk,3));
   end

   statsGwater(jj).tints_p100 = sum([RPixelValues_m100]);
   statsGwater(jj).area_p100 = numel([RPixelValues_m100]);
   statsGwater(jj).volume_p100 = numel([RPixelValues_m100]>0);
   statsGwater(jj).WeightedCentroid_p100(1) = sum([PixelList_m100(:,1)].*...
           double([RPixelValues_m100]))/(sum([RPixelValues_m100]));
   statsGwater(jj).WeightedCentroid_p100(2) = sum([PixelList_m100(:,2)].*...
           double([RPixelValues_m100]))/(sum([RPixelValues_m100]));
   statsGwater(jj).WeightedCentroid_p100(3) =  sum([PixelList_m100(:,3)].*...
           double([RPixelValues_m100]))/(sum([RPixelValues_m100]));
%also add these statxture, 
   statsGwater(jj).statxtureP_100all = statxture2([RPixelValues_m100]);
   posvals2 = [RPixelValues_m100]>0;
   statsGwater(jj).statxtureP_100pos = statxture2([RPixelValues_m100(posvals2)]);
   centG_p100(jj,:) = statsGwater(jj).WeightedCentroid_p100;
   tintsG_p100(jj) = statsGwater(jj).tints_p100;
   volumeG_p100(jj) =statsGwater(jj).volume_p100;
   
      size2= 140;   
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m140 = [x,y,z];
   RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
   for kk=1:numel(PixelList_m140(:,1))
      RPixelValues_m140(kk) = curr1a(PixelList_m140(kk,2),...
          PixelList_m140(kk,1), PixelList_m140(kk,3));
   end

   statsGwater(jj).tints_p140 = sum([RPixelValues_m140]);
      statsGwater(jj).area_p140 = numel([RPixelValues_m140]);

   statsGwater(jj).volume_p140 = numel([RPixelValues_m140]>0);
   statsGwater(jj).WeightedCentroid_p140(1) = sum([PixelList_m140(:,1)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsGwater(jj).WeightedCentroid_p140(2) = sum([PixelList_m140(:,2)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsGwater(jj).WeightedCentroid_p140(3) =  sum([PixelList_m140(:,3)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
%also add these statxture, 
   statsGwater(jj).statxtureP_140all = statxture2([RPixelValues_m140]);
   posvals2 = [RPixelValues_m140]>0;
   statsGwater(jj).statxtureP_140pos = statxture2([RPixelValues_m140(posvals2)]);
   centG_p140(jj,:) = statsGwater(jj).WeightedCentroid_p140;
   tintsG_p140(jj) = statsGwater(jj).tints_p140;
   volumeG_p140(jj) =statsGwater(jj).volume_p140;

    size2= 140;   
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m140 = [x,y,z];
   RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
   for kk=1:numel(PixelList_m140(:,1))
      RPixelValues_m140(kk) = curr1b(PixelList_m140(kk,2),...
          PixelList_m140(kk,1), PixelList_m140(kk,3));
   end

   statsGwater(jj).tints_P140 = sum([RPixelValues_m140]);
   statsGwater(jj).volume_P140 = numel([RPixelValues_m140]>0);
      statsGwater(jj).area_P140 = numel([RPixelValues_m140]);
   statsGwater(jj).WeightedCentroid_P140(1) = sum([PixelList_m140(:,1)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsGwater(jj).WeightedCentroid_P140(2) = sum([PixelList_m140(:,2)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsGwater(jj).WeightedCentroid_P140(3) =  sum([PixelList_m140(:,3)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
%also add these statxture, 
   statsGwater(jj).statxtureP_P140all = statxture2([RPixelValues_m140]);
   posvals2 = [RPixelValues_m140]>0;
   statsGwater(jj).statxtureP_P140pos = statxture2([RPixelValues_m140(posvals2)]);
   centG_P140(jj,:) = statsGwater(jj).WeightedCentroid_P140;
   tintsG_P140(jj) = statsGwater(jj).tints_P140;
   volumeG_P140(jj) =statsGwater(jj).volume_P140;
   
      size2= 180;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m180 = [x,y,z];
   RPixelValues_m180 = zeros(numel(PixelList_m180(:,1)),1);
   for kk=1:numel(PixelList_m180(:,1))
      RPixelValues_m180(kk) = curr1a(PixelList_m180(kk,2),...
          PixelList_m180(kk,1), PixelList_m180(kk,3));
   end

   statsGwater(jj).tints_p180 = sum([RPixelValues_m180]);
   statsGwater(jj).volume_p180 = numel([RPixelValues_m180]>0);
      statsGwater(jj).area_p180 = numel([RPixelValues_m180]);

   statsGwater(jj).WeightedCentroid_p180(1) = sum([PixelList_m180(:,1)].*...
           double([RPixelValues_m180]))/(sum([RPixelValues_m180]));
   statsGwater(jj).WeightedCentroid_p180(2) = sum([PixelList_m180(:,2)].*...
           double([RPixelValues_m180]))/(sum([RPixelValues_m180]));
   statsGwater(jj).WeightedCentroid_p180(3) =  sum([PixelList_m180(:,3)].*...
           double([RPixelValues_m180]))/(sum([RPixelValues_m180]));
%also add these statxture, 
   statsGwater(jj).statxtureP_180all = statxture2([RPixelValues_m180]);
   posvals2 = [RPixelValues_m180]>0;
   statsGwater(jj).statxtureP_180pos = statxture2([RPixelValues_m180(posvals2)]);
   centG_p180(jj,:) = statsGwater(jj).WeightedCentroid_p180;
   tintsG_p180(jj) = statsGwater(jj).tints_p180;
   volumeG_p180(jj) =statsGwater(jj).volume_p180;
    
   
      size2= 220;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m220 = [x,y,z];
   RPixelValues_m220 = zeros(numel(PixelList_m220(:,1)),1);
   for kk=1:numel(PixelList_m220(:,1))
      RPixelValues_m220(kk) = curr1a(PixelList_m220(kk,2),...
          PixelList_m220(kk,1), PixelList_m220(kk,3));
   end

   statsGwater(jj).tints_p220 = sum([RPixelValues_m220]);
   statsGwater(jj).volume_p220 = numel([RPixelValues_m220]>0);
         statsGwater(jj).area_p220 = numel([RPixelValues_m220]);

   statsGwater(jj).WeightedCentroid_p220(1) = sum([PixelList_m220(:,1)].*...
           double([RPixelValues_m220]))/(sum([RPixelValues_m220]));
   statsGwater(jj).WeightedCentroid_p220(2) = sum([PixelList_m220(:,2)].*...
           double([RPixelValues_m220]))/(sum([RPixelValues_m220]));
   statsGwater(jj).WeightedCentroid_p220(3) =  sum([PixelList_m220(:,3)].*...
           double([RPixelValues_m220]))/(sum([RPixelValues_m220]));
%also add these statxture, 
   statsGwater(jj).statxtureP_220all = statxture2([RPixelValues_m220]);
   posvals2 = [RPixelValues_m220]>0;
   statsGwater(jj).statxtureP_220pos = statxture2([RPixelValues_m220(posvals2)]);
   centG_p220(jj,:) = statsGwater(jj).WeightedCentroid_p220;
   tintsG_p220(jj) = statsGwater(jj).tints_p220;
   volumeG_p220(jj) =statsGwater(jj).volume_p220;
    
   
      size2= 260;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
      [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m260 = [x,y,z];
   RPixelValues_m260 = zeros(numel(PixelList_m260(:,1)),1);
   for kk=1:numel(PixelList_m260(:,1))
      RPixelValues_m260(kk) = curr1a(PixelList_m260(kk,2),...
          PixelList_m260(kk,1), PixelList_m260(kk,3));
   end

   statsGwater(jj).tints_p260 = sum([RPixelValues_m260]);
   statsGwater(jj).volume_p260 = numel([RPixelValues_m260]>0);
         statsGwater(jj).area_p260 = numel([RPixelValues_m260]);

   statsGwater(jj).WeightedCentroid_p260(1) = sum([PixelList_m260(:,1)].*...
           double([RPixelValues_m260]))/(sum([RPixelValues_m260]));
   statsGwater(jj).WeightedCentroid_p260(2) = sum([PixelList_m260(:,2)].*...
           double([RPixelValues_m260]))/(sum([RPixelValues_m260]));
   statsGwater(jj).WeightedCentroid_p260(3) =  sum([PixelList_m260(:,3)].*...
           double([RPixelValues_m260]))/(sum([RPixelValues_m260]));
%also add these statxture, 
   statsGwater(jj).statxtureP_260all = statxture2([RPixelValues_m260]);
   posvals2 = [RPixelValues_m260]>0;
   statsGwater(jj).statxtureP_260pos = statxture2([RPixelValues_m260(posvals2)]);
   centG_p260(jj,:) = statsGwater(jj).WeightedCentroid_p260;
   tintsG_p260(jj) = statsGwater(jj).tints_p260;
   volumeG_p260(jj) =statsGwater(jj).volume_p260;
    
   
      size2= 300;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m300 = [x,y,z];
   RPixelValues_m300 = zeros(numel(PixelList_m300(:,1)),1);
   for kk=1:numel(PixelList_m300(:,1))
      RPixelValues_m300(kk) = curr1a(PixelList_m300(kk,2),...
          PixelList_m300(kk,1), PixelList_m300(kk,3));
   end

   statsGwater(jj).tints_p300 = sum([RPixelValues_m300]);
   statsGwater(jj).volume_p300 = numel([RPixelValues_m300]>0);
         statsGwater(jj).area_p300 = numel([RPixelValues_m300]);

   statsGwater(jj).WeightedCentroid_p300(1) = sum([PixelList_m300(:,1)].*...
           double([RPixelValues_m300]))/(sum([RPixelValues_m300]));
   statsGwater(jj).WeightedCentroid_p300(2) = sum([PixelList_m300(:,2)].*...
           double([RPixelValues_m300]))/(sum([RPixelValues_m300]));
   statsGwater(jj).WeightedCentroid_p300(3) =  sum([PixelList_m300(:,3)].*...
           double([RPixelValues_m300]))/(sum([RPixelValues_m300]));
%also add these statxture, 
   statsGwater(jj).statxtureP_300all = statxture2([RPixelValues_m300]);
   posvals2 = [RPixelValues_m300]>0;
   statsGwater(jj).statxtureP_300pos = statxture2([RPixelValues_m300(posvals2)]);
   centG_p300(jj,:) = statsGwater(jj).WeightedCentroid_p300;
   tintsG_p300(jj) = statsGwater(jj).tints_p300;
   volumeG_p300(jj) =statsGwater(jj).volume_p300;
    
      size2= 400;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m400 = [x,y,z];
   RPixelValues_m400 = zeros(numel(PixelList_m400(:,1)),1);
   for kk=1:numel(PixelList_m400(:,1))
      RPixelValues_m400(kk) = curr1a(PixelList_m400(kk,2),...
          PixelList_m400(kk,1), PixelList_m400(kk,3));
   end

   statsGwater(jj).tints_p400 = sum([RPixelValues_m400]);
   statsGwater(jj).volume_p400 = numel([RPixelValues_m400]>0);
         statsGwater(jj).area_p400 = numel([RPixelValues_m400]);

   statsGwater(jj).WeightedCentroid_p400(1) = sum([PixelList_m400(:,1)].*...
           double([RPixelValues_m400]))/(sum([RPixelValues_m400]));
   statsGwater(jj).WeightedCentroid_p400(2) = sum([PixelList_m400(:,2)].*...
           double([RPixelValues_m400]))/(sum([RPixelValues_m400]));
   statsGwater(jj).WeightedCentroid_p400(3) =  sum([PixelList_m400(:,3)].*...
           double([RPixelValues_m400]))/(sum([RPixelValues_m400]));
%also add these statxture, 
   statsGwater(jj).statxtureP_400all = statxture2([RPixelValues_m400]);
   posvals2 = [RPixelValues_m400]>0;
   statsGwater(jj).statxtureP_400pos = statxture2([RPixelValues_m400(posvals2)]);
   centG_p400(jj,:) = statsGwater(jj).WeightedCentroid_p400;
   tintsG_p400(jj) = statsGwater(jj).tints_p400;
   volumeG_p400(jj) =statsGwater(jj).volume_p400;
    
   
        size2= 500;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m500 = [x,y,z];
   RPixelValues_m500 = zeros(numel(PixelList_m500(:,1)),1);
   for kk=1:numel(PixelList_m500(:,1))
      RPixelValues_m500(kk) = curr1a(PixelList_m500(kk,2),...
          PixelList_m500(kk,1), PixelList_m500(kk,3));
   end

   statsGwater(jj).tints_p500 = sum([RPixelValues_m500]);
   statsGwater(jj).volume_p500 = numel([RPixelValues_m500]>0);
         statsGwater(jj).area_p500 = numel([RPixelValues_m500]);

   statsGwater(jj).WeightedCentroid_p500(1) = sum([PixelList_m500(:,1)].*...
           double([RPixelValues_m500]))/(sum([RPixelValues_m500]));
   statsGwater(jj).WeightedCentroid_p500(2) = sum([PixelList_m500(:,2)].*...
           double([RPixelValues_m500]))/(sum([RPixelValues_m500]));
   statsGwater(jj).WeightedCentroid_p500(3) =  sum([PixelList_m500(:,3)].*...
           double([RPixelValues_m500]))/(sum([RPixelValues_m500]));
%also add these statxture, 
   statsGwater(jj).statxtureP_500all = statxture2([RPixelValues_m500]);
   posvals2 = [RPixelValues_m500]>0;
   statsGwater(jj).statxtureP_500pos = statxture2([RPixelValues_m500(posvals2)]);
   centG_p500(jj,:) = statsGwater(jj).WeightedCentroid_p500;
   tintsG_p500(jj) = statsGwater(jj).tints_p500;
   volumeG_p500(jj) =statsGwater(jj).volume_p500;
    
   
   statsGwater(jj).statxture_all = statxture2([statsGwater(jj).PixelValues]);
   posvals = [statsGwater(jj).PixelValues]>0;
   statsGwater(jj).statxture_pos = statxture2([statsGwater(jj).PixelValues(posvals)]);
   statxture_all(jj,:) = statsGwater(jj).statxture_all;
   statxture_pos(jj,:) = statsGwater(jj).statxture_pos;
end
%
save([path 'add_to_statsGw10_edges.mat'],'centG_p20','tintsG_p20',...)
    'volumeG_p20','centG_p60','tintsG_p60',...)
    'volumeG_p60','centG_p100','tintsG_p100',...)
    'volumeG_p100','centG_p140','tintsG_p140',...)
    'volumeG_p140','centG_p180','tintsG_p180',...)
    'volumeG_p180','centG_p220','tintsG_p220',...)
    'volumeG_p220','centG_p260','tintsG_p260',...)
    'volumeG_p260','centG_p300','tintsG_p300',...)
    'volumeG_p300','centG_p400','tintsG_p400',...)
    'volumeG_p400','centG_p500','tintsG_p500',...)
    'volumeG_p500','centG_P140','tintsG_P140',...)
    'volumeG_P140','statxture_all','statxture_pos','-v7.3')
save([path 'statsG2w10_edges_plus.mat'],'statsGwater','-v7.3')