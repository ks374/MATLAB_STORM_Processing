clear all
clc
%
base_path = 'Y:\Chenghang\04_4_Color\chenghaz_001_Sample_01_A\';
exp_folder = [base_path 'analysis\'];
path = [exp_folder  'elastic_align\'];

mergedpath = [path 'storm_merged\'];
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);


convpath = [path 'conv_merged\'];
convfiles = [dir([convpath '*.png']) dir([convpath '*.tif'])];

mask_folder = [exp_folder 'Result\1_Soma\'];
maskfiles = [dir([mask_folder '*.tif']) dir([mask_folder '*.png'])];

%set voxel size in nm
voxel=[15.5, 15.5, 70];
outpath = [exp_folder 'Result\'];

%
clear BG
disp('allocating arrays')
BG = zeros(info.Height, info.Width, num_images,'uint8');
disp('loading data')
parfor k = 1:num_images
    A = imread([convpath convfiles(k,1).name]);
    BG(:,:,k) = A(:,:,1);
end

% bad_sec = [4,17];
% BG(:,:,bad_sec) = [];
num_images_2 = size(BG,3);

parfor k =1:num_images_2
    M = imread([mask_folder maskfiles(k,1).name]);
    Mask(:,:,k) = M(:,:);
end
%

Mask = logical(Mask);
for i = 1:num_images_2
    I(:,:,i) = BG(:,:,i).*uint8(1.-Mask(:,:,i));
%     Masked_out(:,:,i) = BG(:,:,i).*uint8(Mask(:,:,i));
end
%
% parfor i = 1:num_images
%     imwrite(Masked_out(:,:,i),[outpath 'Maked_outR_' sprintf('%03d',i) '.tif']);
% end
%
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
%frame = zeros(info.Height, info.Width, num_images,'uint8');
DYs = false(size(BG));
disp('thresholding YFP')
parfor k = 1:num_images_2
    DYs(:,:,k) = im2bw(I(:,:,k), leveluse);
end
%
%     Result(:,:,1) = I(:,:,1).*uint8(DYs(:,:,1));
%     imwrite(Result(:,:,1),[outpath 'Maked_' sprintf('%03d',1) '.tif']);

%

S = zeros(info.Height, info.Width, num_images,'uint8');
disp('loading data')
parfor k = 1:num_images
    A = imread([mergedpath mergedfiles(k,1).name]);
    S(:,:,k) = A(:,:,1);
end
% S(:,:,bad_sec) = [];

%
for i = 1:num_images_2
    S1(:,:,i) = S(:,:,i).*uint8(DYs(:,:,i));
end
%
% parfor i = 1:num_images
%     imwrite(S1(:,:,i),[outpath 'Maked_G_storm_' sprintf('%03d',i) '.tif']);
% end

%
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

threshfactorg = double(multithresh(hy1dist,2));
t_use = threshfactorg(2)/256

%
disp('making CG')
CG = false(size(S2));
parfor k=1:size(S2,3)
    CG(:,:,k) = im2bw(S2(:,:,k), t_use);
end
%
% parfor i = 1:num_images
%     imwrite(uint8(CG(:,:,i)),[outpath 'mask2_' sprintf('%03d',i) '.tif']);
% end

%
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
areacutoff =4 ;
statsGgauss = statsGgauss([statsG.Area]>areacutoff);
statsG = statsG([statsG.Area]>areacutoff);
%
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

% save([outpath 'sizeshapematG_and_cent_water10_area_cutoff.mat'],'centGw','sizeshape_mat')
% save([outpath 'statsGwater10_area_cutoff.mat'],'statsGwater','-v7.3')
%
No = [];
for i = 1:numel(statsGwater)
    minpix = min(statsGwater(i).PixelIdxList);
    maxpix = max(statsGwater(i).PixelIdxList);
    if Mask(minpix) == 1 | Mask(maxpix) == 1
        No = cat(1,No,i);
    end
end
statsGwater(No) = [];
centGw(No,:) = [];
sizeshape_mat(No,:) = [];

%
centG = centGw;
sizeshape_matG = sizeshape_mat;

for i=1:numel(centG(:,1))
   centG(i,:) = centG(i,:).*voxel ;
end
rcentG = centG;
%%
i=4;  % i = range(8)
sizeshapuse = sizeshape_mat;
val1w = (sizeshapuse(:,10+i)+1)./(sizeshapuse(:,2+i)+1);
val2w = log10((sizeshapuse(:,10+i)+1)*0.0158*0.0158*0.07);
val2w = log((sizeshapuse(:,10+i)+1));
numel(find(val1w<0.99))

Xn=80; Yn=60;
Xrange=[min(val1w) max(val1w)]; Yrange=[min(val2w) max(val2w)];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)'; 
%
figure;
H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
close

cutoffg = 60;
H2 = H;
H2(H2>cutoffg)=cutoffg;

figure;
pcolor(X,Y,H2)
%
savefig([outpath 'Gsizeshape_heatmap.fig'])
%%

figure;
pcolor(X,Y,H2)
% manually draw polygon on figure
currpoly=impoly
% return polygon coordinates
synapse_regiong=currpoly.getPosition
% save figure of selected polygon region
savefig([outpath 'synapse_selection_poly.fig'])

%Return centroid and stats lists for all clusters in selected polygon area
centGa2s=centG(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

centGa2ns=centG(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

rcentGa2s=rcentG(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

rcentGa2ns=rcentG(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

% return size/shape array for selected clusters
sizeshape_matGa2s=sizeshape_matG(find(inpolygon(val1w,val2w,...
    synapse_regiong(:,1),synapse_regiong(:,2))),:);

sizeshape_matGa2ns=sizeshape_matG(find(~inpolygon(val1w,val2w,...
    synapse_regiong(:,1),synapse_regiong(:,2))),:);
 
% return detailed pixel stats lists for selected clusters
statsGa2s=statsGwater(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

statsGa2ns=statsGwater(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

size(centGa2s,1)

%
Real_non_select_area_thre = 4;
centGa2ns = centGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
rcentGa2ns = rcentGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
statsGa2ns = statsGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
sizeshape_matGa2ns = sizeshape_matGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);

%save([base_folder 'statslistG2sw10.mat'],'centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statslistG2nsw10.mat'],'centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
save([outpath 'statslistR2sw10.mat'],'centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statslistR2nsw10.mat'],'centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')

%
save([outpath 'statsR2sw10.mat'],'statsGa2s','centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statsR2nsw10.mat'],'statsGa2ns','centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')
%save([base_folder 'statsG2sw10.mat'],'statsGb2s','centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statsG2nsw10.mat'],'statsGb2ns','centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
