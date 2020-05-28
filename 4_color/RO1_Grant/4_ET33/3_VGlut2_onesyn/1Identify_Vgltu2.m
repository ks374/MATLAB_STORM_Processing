clear all
clc
%%
base_path = 'Z:\Chenghang\chenghaz_008_Sample_Epi_Ipsi_A\';
exp_folder = [base_path 'analysis\'];
path = [exp_folder  'elastic_align\'];

mergedpath = [path 'storm_merged\'];
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);


convpath = [path 'conv_merged\'];
convfiles = [dir([convpath '*.png']) dir([convpath '*.tif'])];

mask_folder = [path 'Result/1_Soma/'];
maskfiles = [dir([mask_folder '*.tif']) dir([mask_folder '*.png'])];

%set voxel size in nm
voxel=[15.5, 15.5, 70];
outpath = 'Z:\Chenghang\chenghaz_008_Sample_Epi_Ipsi_A\analysis\elastic_align\Result\VGlut2\';

%%
clear BG
disp('allocating arrays')
BG = zeros(info.Height, info.Width, num_images,'uint8');
disp('loading data')
parfor k = 1:num_images
    A = imread([convpath convfiles(k,1).name]);
    BG(:,:,k) = A(:,:,3);
end
% 
% bad_sec = [2];
% BG(:,:,bad_sec) = [];
num_images_2 = size(BG,3);

parfor k=1:num_images_2
    M = imread([mask_folder maskfiles(k,1).name]);
    Mask(:,:,k) = M(:,:);
end

%

Mask = logical(Mask);
for i = 1:num_images_2
    I(:,:,i) = BG(:,:,i).*uint8(1.-Mask(:,:,i));
    Masked_out(:,:,i) = BG(:,:,i).*uint8(Mask(:,:,i));
end
%%

parfor i = 1:num_images
    imwrite(Masked_out(:,:,i),[outpath 'Maked_outR_' sprintf('%03d',i) '.tif']);
end
%%
clear hy1
parfor j=1:num_images_2
    disp(j)
     A1 = I(:,:,j); 
     A1a = A1(find(A1));
     [hy, hx] = hist(A1a,0:1:255);
     hy1(j,:) = hy;
end
%
hx2 = 0:255;
hy1dist = [];
mhy1 = mean(hy1)./10;
for i=1:numel(mhy1)
nums = ones(1,round(mhy1(i))).*hx2(i);       
hy1dist = cat(2,hy1dist,nums);
end

levelYs = multithresh(hy1dist,2);
leveluse = double(levelYs(1))/256

%
DYs = false(size(BG));
disp('thresholding YFP')
parfor k = 1:num_images_2
    DYs(:,:,k) = im2bw(I(:,:,k), leveluse);
end
%%
% parfor i = 1:num_images
%     Result(:,:,i) = I(:,:,i).*uint8(DYs(:,:,i));
% end
% 
% parfor i = 1:num_images
%     imwrite(Result(:,:,i),[outpath 'Maked_' sprintf('%03d',i) '.tif']);
% end

%%

S = zeros(info.Height, info.Width, num_images,'uint8');
disp('loading data')
parfor k = 1:num_images
    A = imread([mergedpath mergedfiles(k,1).name]);
    S(:,:,k) = A(:,:,3);
end
% S(:,:,bad_sec) = [];
%
for i = 1:num_images_2
    S1(:,:,i) = S(:,:,i).*uint8(DYs(:,:,i));
end
%%
parfor i = 1:10
    imwrite(S1(:,:,i),[outpath 'Maked_S_' sprintf('%03d',i) '.tif']);
end

%%
gauss  = 5;
gausspix = (gauss);
sigmaZ = gausspix* voxel(1)/voxel(3);
sizZ= gausspix*2*sigmaZ* voxel(1)/voxel(3);
xZ=-ceil(sizZ/2):ceil(sizZ/2);
H1 = exp(-(xZ.^2/(2*sigmaZ^2)));
H1 = H1/sum(H1(:));
Hz=reshape(H1,[1 1 length(H1)]);
%create blurred gephyrin for region
%
S2 = zeros(size(BG),'uint8');
parfor k=1:num_images_2
    disp(k)
    S2(:,:,k) = (imfilter(S1(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
end

S2 = imfilter(S2,Hz,'same','replicate');
%
clear hy1
parfor j=1:num_images_2
    disp(j)
     A1 = S2(:,:,j); 
     A1a = A1(find(A1));
     [hy, hx] = hist(A1a,0:1:255);
     hy1(j,:) = hy;
end
%
hx2 = 0:255;
hy1dist = [];
mhy1 = mean(hy1)./10;
for i=1:numel(mhy1)
nums = ones(1,round(mhy1(i))).*hx2(i);       
hy1dist = cat(2,hy1dist,nums);
end
%
% figure;
% hist(hy1dist,255)

threshfactorg = double(multithresh(hy1dist,2))
t_use = threshfactorg(2)/256

%
disp('making CG')
CG = false(size(S2));
parfor k=1:size(S2,3)
    CG(:,:,k) = im2bw(S2(:,:,k), t_use);
end
%%
% parfor i = 1:num_images
%     imwrite(uint8(CG(:,:,i)),[outpath 'mask2_' sprintf('%03d',i) '.tif']);
% end

%%
%clear bg2
disp('making CCG')
CCG = bwconncomp(CG,26); 
%clear CG
disp('making statsG')
%
statsG = regionprops(CCG,S1,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
statsGgauss = regionprops(CCG,S2,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');

%
statsG_backup = statsG;
statsGgauss_backup = statsGgauss;

%
areacutoff = 4;
statsGgauss = statsGgauss(log(1+[statsG.Area])>areacutoff);
statsG = statsG(log(1+[statsG.Area])>areacutoff);
%%

% sizeBG = size(BG);
% clear CG CCG
% clear statsG2a
% 
% statsG2a(numel(statsGgauss)).stats = [];
% clear minpix maxpix curr2 curr2b I2 I3 L L2 curr3 x y z
% 
% parfor jj=1:numel(statsGgauss)
%    disp(jj)
% minpix = min(statsGgauss(jj).PixelList);
% maxpix = max(statsGgauss(jj).PixelList);
% curr2 = zeros((maxpix(2)-minpix(2)+1+4),(maxpix(1)-minpix(1)+1+4),maxpix(3)- ...
% minpix(3)+1+2,'uint8');
% curr2b = zeros((maxpix(2)-minpix(2)+1+4),(maxpix(1)-minpix(1)+1+4),maxpix(3)- ...
% minpix(3)+1+2,'uint8');
% for j=1: numel(statsGgauss(jj).PixelList(:,1))
%    curr2(statsG(jj).PixelList(j,2)-minpix(2)+1+2,statsG(jj).PixelList(j,1)-minpix(1)+1+2, ...
%        statsG(jj).PixelList(j,3)-minpix(3)+1+1)=statsG(jj).PixelValues(j,1); 
%    curr2b(statsGgauss(jj).PixelList(j,2)-minpix(2)+1+2,statsGgauss(jj).PixelList(j,1)-minpix(1)+1+2, ...
%        statsGgauss(jj).PixelList(j,3)-minpix(3)+1+1)=statsGgauss(jj).PixelValues(j,1); 
% end
% 
% % parfor i = 1:num_images
% %     imwrite(S1(:,:,i),[outpath 'Maked_S_' sprintf('%03d',i) '.tif']);
% % end
% 
% 
% I2 = imcomplement(curr2b);
% I3 = imhmin(I2,10);
% L = watershed(I3);
% for i=1:max(L(:))
% L2 = L;
% curr3 = curr2b;
% curr3(L~=i) = 0;
% L2(L~=i) = 0;
% 
% [y, x, z] = ind2sub(size(curr3),find(curr3));
% statsG2a(jj).stats(i).PixelList = [x y z];
% for kk=1:numel(statsG2a(jj).stats(i).PixelList(:,1))
% statsG2a(jj).stats(i).PixelValues(kk,1) = ...
%     curr2(statsG2a(jj).stats(i).PixelList(kk,2), ...
%     statsG2a(jj).stats(i).PixelList(kk,1),...
%     statsG2a(jj).stats(i).PixelList(kk,3));
% end
% 
% statsG2a(jj).stats(i).PixelList(:,1)= ...
%     statsG2a(jj).stats(i).PixelList(:,1) + minpix(1) - 1-2;
% statsG2a(jj).stats(i).PixelList(:,2)= ...
%     statsG2a(jj).stats(i).PixelList(:,2) + minpix(2) - 1-2;
% statsG2a(jj).stats(i).PixelList(:,3)= ...
%     statsG2a(jj).stats(i).PixelList(:,3) + minpix(3) - 1-1;
% statsG2a(jj).stats(i).Area = numel(statsG2a(jj).stats(i).PixelValues);
% statsG2a(jj).stats(i).TintsG = sum(statsG2a(jj).stats(i).PixelValues);
% statsG2a(jj).stats(i).PixelIdxList = ...
%     sub2ind(sizeBG,statsG2a(jj).stats(i).PixelList(:,2),...
%     statsG2a(jj).stats(i).PixelList(:,1),...
%     statsG2a(jj).stats(i).PixelList(:,3));
% end
% end
%%
%%
% parfor i = 1:size(curr2,3)
%      imwrite(double(W(:,:,i)),[outpath 'W_' sprintf('%03d',i) '.tif']);
% end%%

%%
% clear templist
% templist = []; statsG2a2 = []; cnt=0;
% for i=1:numel(statsG2a)
% templist = cat(1,templist, statsG2a(i).stats');
% cnt = cnt+1;
%      if cnt==2000
% statsG2a2 = cat(1,statsG2a2, templist);
% templist = [];
% cnt = 0;
%      end
% end
% statsG2a2 = cat(1,statsG2a2, templist);
% statsGwater = statsG2a2;
%%

statsGwater = statsG;

for i=1:numel(statsGwater)
    disp(i)
    statsGwater(i,1).PixelValues = S1(statsGwater(i,1).PixelIdxList);
    statsGwater(i,1).PixelValues2 = S2(statsGwater(i,1).PixelIdxList);
statsGwater(i,1).Volume1_0 = numel(find(statsGwater(i,1).PixelValues>0));
statsGwater(i,1).Volume2_0 = numel(find(statsGwater(i,1).PixelValues2>0));
statsGwater(i,1).Volume2_t2a = numel(find(statsGwater(i,1).PixelValues2>(threshfactorg(1))));
statsGwater(i,1).Volume2_t2b = numel(find(statsGwater(i,1).PixelValues2>(1.1*threshfactorg(1))));
statsGwater(i,1).Volume2_t2c = numel(find(statsGwater(i,1).PixelValues2>(1.2*threshfactorg(1))));
statsGwater(i,1).Volume2_t2d = numel(find(statsGwater(i,1).PixelValues2>(1.3*threshfactorg(1))));
statsGwater(i,1).Volume2_t2e = numel(find(statsGwater(i,1).PixelValues2>(1.4*threshfactorg(1))));
statsGwater(i,1).Volume2_t2f = numel(find(statsGwater(i,1).PixelValues2>(1.5*threshfactorg(1))));
statsGwater(i,1).Volume2_t2g = numel(find(statsGwater(i,1).PixelValues2>(1.7*threshfactorg(1))));
statsGwater(i,1).Volume2_t2h = numel(find(statsGwater(i,1).PixelValues2>(2.0*threshfactorg(1))));
statsGwater(i,1).Volume1_t2a0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2b0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.1*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2c0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.2*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2d0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.3*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2e0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.4*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2f0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.5*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2g0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.7*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2h0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    2.0*threshfactorg(1))>0));

statsGwater(i,1).Area = numel(statsGwater(i,1).PixelValues);
statsGwater(i,1).TintsG = sum(statsGwater(i,1).PixelValues);

statsGwater(i).WeightedCentroid(1) = ...
        sum([statsGwater(i).PixelList(:,1)].*...
        double([statsGwater(i).PixelValues]))/...
        (sum([statsGwater(i).PixelValues]));
statsGwater(i).WeightedCentroid(2) = ...
        sum([statsGwater(i).PixelList(:,2)].*...
        double([statsGwater(i).PixelValues]))/...
        (sum([statsGwater(i).PixelValues]));
statsGwater(i).WeightedCentroid(3) = ...
        sum([statsGwater(i).PixelList(:,3)].*...
        double([statsGwater(i).PixelValues]))/...
        (sum([statsGwater(i).PixelValues]));
end

numel(find([statsGwater.TintsG]>0))
statsGwater = statsGwater([statsGwater.TintsG]>0);
%
sizeshape_mat = cat(2,[statsGwater.Volume1_0]',...
    [statsGwater.Volume2_0]',[statsGwater.Volume2_t2a]',[statsGwater.Volume2_t2b]',...
    [statsGwater.Volume2_t2c]',[statsGwater.Volume2_t2d]',[statsGwater.Volume2_t2e]',...
    [statsGwater.Volume2_t2f]',[statsGwater.Volume2_t2g]',[statsGwater.Volume2_t2h]',...
    [statsGwater.Volume1_t2a0]', [statsGwater.Volume1_t2b0]',...
    [statsGwater.Volume1_t2c0]',[statsGwater.Volume1_t2d0]',...
    [statsGwater.Volume1_t2e0]',[statsGwater.Volume1_t2f0]',...
    [statsGwater.Volume1_t2g0]', [statsGwater.Volume1_t2h0]',...
    [statsGwater.Area]', [statsGwater.TintsG]');

centGw = zeros(numel(statsGwater),3);
dsvoxel = 155/0.3;
%dsvoxel = 158;
voxel = [15.5 15.5 70];
%
for jj =1:numel(statsGwater)   
    centGw(jj,:) = statsGwater(jj,1).WeightedCentroid;
end

%

save([outpath 'sizeshapematV_and_cent_water10_area_cutoff.mat'],'centGw','sizeshape_mat')
save([outpath 'statsVwater10_area_cutoff.mat'],'statsGwater','-v7.3')
%%
clear