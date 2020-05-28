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
load([path 'statsG2sw10.mat']);
%%
clear BP BP2
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsGa2s)
    disp(int2str(i))
    BP(statsGa2s(i).PixelIdxList)=statsGa2s(i).PixelValues;  
end
%BP2 = BP;
%for i = 1:numel(statsGa2ns)
%    disp(int2str(i))
%    BP2(statsGa2ns(i).PixelIdxList)=statsGa2ns(i).PixelValues;  
%end
%%
load([path 'statsP2sw10.mat']);
%%
statsGwater = statsPa2s;
%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsGwater):-1:1
    disp(i)
statsGwater(i).tints_p100 = [];
statsGwater(i).volume_p100 = [];
statsGwater(i).WeightedCentroid_p100 = [];
statsGwater(i).statxtureP_100all = [];
statsGwater(i).statxtureP_100pos = [];
centG_p100(i,1:3) = 0;
tintsG_p100(i) = 0;
volumeG_p100(i) =0;

statsGwater(i).tints_p140 = [];
statsGwater(i).volume_p140 = [];
statsGwater(i).WeightedCentroid_p140 = [];
statsGwater(i).statxtureP_140all = [];
statsGwater(i).statxtureP_140pos = [];
centG_p140(i,1:3) = 0;
tintsG_p140(i) = 0;
volumeG_p140(i) =0;

statsGwater(i).tints_p70 = [];
statsGwater(i).volume_p70 = [];
statsGwater(i).WeightedCentroid_p70 = [];
statsGwater(i).statxtureP_70all = [];
statsGwater(i).statxtureP_70pos = [];
centG_p70(i,1:3) = 0;
tintsG_p70(i) = 0;
volumeG_p70(i) =0;

statsGwater(i).statxture_all = [];
statsGwater(i).statxture_pos = [];
statxture_all(i,1:6) = 0;
statxture_pos(i,1:6) = 0;

%try preproccessing the BP varialbe in a for loop
   minpix = min(statsGwater(i).PixelList);  maxpix = max(statsGwater(i).PixelList);
   min1 = minpix(1)-10; min2 = minpix(2)-10; min3 = minpix(3)-3;
   max1 = maxpix(1)+10; max2 = maxpix(2)+10; max3 = maxpix(3)+3;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end      
   BPp(i).mat = BP(min2:max2,min1:max1,min3:max3); 
end
%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
%
parfor jj=1:numel(statsGwater)
   %
   disp(jj)
   minpix = min(statsGwater(jj).PixelList);  maxpix = max(statsGwater(jj).PixelList);
   min1 = minpix(1)-10; min2 = minpix(2)-10; min3 = minpix(3)-3;
   max1 = maxpix(1)+10; max2 = maxpix(2)+10; max3 = maxpix(3)+3;
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
   Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);

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
    
   
      size2= 70;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
     [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m70 = [x,y,z];
   RPixelValues_m70 = zeros(numel(PixelList_m70(:,1)),1);
   for kk=1:numel(PixelList_m70(:,1))
      RPixelValues_m70(kk) = curr1a(PixelList_m70(kk,2),...
          PixelList_m70(kk,1), PixelList_m70(kk,3));
   end

   statsGwater(jj).tints_p70 = sum([RPixelValues_m70]);
   statsGwater(jj).volume_p70 = numel([RPixelValues_m70]>0);
   statsGwater(jj).WeightedCentroid_p70(1) = sum([PixelList_m70(:,1)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   statsGwater(jj).WeightedCentroid_p70(2) = sum([PixelList_m70(:,2)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   statsGwater(jj).WeightedCentroid_p70(3) =  sum([PixelList_m70(:,3)].*...
           double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
%also add these statxture, 
   statsGwater(jj).statxtureP_70all = statxture2([RPixelValues_m70]);
   posvals2 = [RPixelValues_m70]>0;
   statsGwater(jj).statxtureP_70pos = statxture2([RPixelValues_m70(posvals2)]);
   centG_p70(jj,:) = statsGwater(jj).WeightedCentroid_p70;
   tintsG_p70(jj) = statsGwater(jj).tints_p70;
   volumeG_p70(jj) =statsGwater(jj).volume_p70;
    
   
   
   statsGwater(jj).statxture_all = statxture2([statsGwater(jj).PixelValues]);
   posvals = [statsGwater(jj).PixelValues]>0;
   statsGwater(jj).statxture_pos = statxture2([statsGwater(jj).PixelValues(posvals)]);
   statxture_all(jj,:) = statsGwater(jj).statxture_all;
   statxture_pos(jj,:) = statsGwater(jj).statxture_pos;
end
%
statsPwater = statsGwater;
save([path 'add_to_statsP2w10.mat'],'centG_p100','tintsG_p100',...)
    'volumeG_p100','centG_p140','tintsG_p140',...)
    'volumeG_p140','centG_p70','tintsG_p70',...)
    'volumeG_p70','statxture_all','statxture_pos','-v7.3')
save([path 'statsP2w10_plus.mat'],'statsGwater','-v7.3')

%%
figure;
hist(log([statsGwater.tints_p260]),200)