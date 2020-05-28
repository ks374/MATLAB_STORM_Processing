clear all
%%
exp_folder = 'Z:\Chenghang\chenghaz_008_Sample_Epi_Ipsi_A/analysis/';
path = [exp_folder  'elastic_align/'];

convpath = [path 'storm_merged/'];
convpathds = [path 'storm_merged_ds/'];
convfiles = [dir([convpath '*.tif']) dir([convpath '*.png'])];
convfilesds = [dir([convpathds '*.tif']) dir([convpathds '*.png'])];
num_images = numel(convfiles);
infos = imfinfo([convpath convfiles(1,1).name]);

outpath = [path 'Result/'];

downsample = 10;
%set voxel size in nm
voxels=[(15.5), (15.5), 70];

BYs = zeros(ceil((infos(1,1).Height)/1),ceil((infos(1,1).Width)/1),num_images,'uint8');

parfor k = 1:num_images
    A = imread([convpath convfiles(k,1).name]);
    BYs(:,:,k) = A(:,:,1);
end
%%
% bad_sec = [];
% BYs(:,:,bad_sec) = [];
num_images_2 = size(BYs,3);
%%
gauss  = 5;
gausspix = (gauss);

bg2 = zeros(size(BYs),'uint8');
parfor k=1:num_images_2
    disp(k)
    bg2(:,:,k) = (imfilter(BYs(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
end


clear hy1
parfor j=1:num_images_2
    disp(j)
     A1 = bg2(:,:,j); 
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
t_use = threshfactorg(1)/256;
% t_use = threshfactorg(2)/256;
%%
disp('making CG')
CG = false(size(bg2));
parfor k=1:size(bg2,3)
    CG(:,:,k) = im2bw(bg2(:,:,k), t_use);
end
%
% imwrite(double(CG(:,:,1)),[outpath 'mask' '.tif']);
%
disp('making CCG')
CCG = bwconncomp(CG,26); 
%clear CG
disp('making statsG')

statsG = regionprops(CCG,BYs,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
statsGgauss = regionprops(CCG,bg2,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');

%
new_G = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
% eroded_G = new_G;
statsG_temp = statsGgauss(find(log(1+[statsGgauss.Area])> 10.5));
numel(statsG_temp)
%%
for i = 1:size(statsG_temp,1)
% for i = 1:10
    disp(i)
    PixelList = statsG_temp(i).PixelList;
    PixelValues = statsG_temp(i).PixelValues;
    for j = 1:size(PixelList,1)
        y = PixelList(j,1);
        x = PixelList(j,2);
        z = PixelList(j,3);
        new_G(x,y,z) = PixelValues(j);
    end
end
%%
% imagesc(new_G(:,:,5))
%%
gauss  = 90;
gausspix = (gauss);

bg3 = zeros(size(BYs),'uint8');
parfor k=1:num_images_2
    disp(k)
    bg3(:,:,k) = (imfilter(new_G(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
end


for_fill = logical(bg3);
for i = 1:num_images_2
    fill(:,:,i) = imfill(for_fill(:,:,i),4,'holes');
end

% imshow(fill(:,:,1))
%
gauss  = 60;
gausspix = (gauss);

bg4 = zeros(size(BYs),'uint8');
parfor k=1:num_images_2
    disp(k)
    I = edge(fill(:,:,k));
    I = uint8(I)*255;
    bg4(:,:,k) = imfilter(I,fspecial('gaussian',gausspix*2,gausspix),'same','replicate');
end
bg4 = logical(bg4);
bg4_re = 1.-bg4;
F = fill.*bg4_re;

%%
for i = 1:num_images_2
    imwrite(double(F(:,:,i)),[outpath 'F_' sprintf('%03d',i) '.tif']);
end