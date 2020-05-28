clear all
clc
%
exp_folder = 'Y:\Chenghang\04_4_Color\chenghaz_005_Sample_P8_C_A\analysis\';
Box_thick = 20;
Dis_thick = 60;
%
path = [exp_folder  'elastic_align/'];

stormpath = [path 'storm_merged/'];
stormfiles = [dir([stormpath '*.tif']) dir([stormpath '*.png'])];
num_images = numel(stormfiles);
infos = imfinfo([stormpath stormfiles(1,1).name]);

convpath = [path 'conv_merged\'];
convfiles = [dir([convpath '*.tif']) dir([convpath '*.png'])];

mkdir([exp_folder 'Result']);
mkdir([exp_folder 'Result\1_soma']);
mkdir([exp_folder 'Result\0_BD'])
outpath = [exp_folder 'Result\1_Soma\'];
voxels=[(15.5), (15.5), 70];

BYs = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
BD = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
Area = zeros(num_images,1);
%
parfor k = 1:num_images
    disp(k)
    A = imread([convpath convfiles(k,1).name]);
    C = A(:,:,2) * 255;
    gausspix = 3;
    B = imfilter(C,fspecial('gaussian',gausspix,gausspix),'same','replicate');
    B = B * 255 - C;
    gausspix = 2*Box_thick + 1;
    D = imfilter(B,fspecial('gaussian',gausspix,gausspix),'same','replicate');
    D = D * 255;
    BD(:,:,k) = D;
end

%
parfor k = 1:num_images
    A = imread([stormpath stormfiles(k,1).name]);
    BYs(:,:,k) = A(:,:,2);
end

% bad_sec = [2];
% BYs(:,:,bad_sec) = [];
num_images_2 = size(BYs,3);
%
gauss  = 10;
gausspix = (gauss);

bg2 = zeros(size(BYs),'uint8');
parfor k=1:num_images_2
    disp(k)
    bg2(:,:,k) = (imfilter(BYs(:,:,k),fspecial('gaussian',gausspix,gausspix),'same','replicate'));
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

threshfactorg = double(multithresh(hy1dist,3));
t_use = threshfactorg(1)/256
% t_use = threshfactorg(2)/256;
%
disp('making CG')
CG = false(size(bg2));
parfor k=1:size(bg2,3)
    CG(:,:,k) = im2bw(bg2(:,:,k), t_use);
end
% imwrite(double(CG(:,:,1)),[outpath 'mask' '.tif']);
%
%
disp('making CCG')
CCG = bwconncomp(CG,26); 
%clear CG
disp('making statsG')

statsG = regionprops(CCG,BYs,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
statsGgauss = regionprops(CCG,bg2,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');

%
new_G = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
%
% eroded_G = new_G;
statsG_temp = statsGgauss(find(log(1+[statsGgauss.Area])> 11.5));
numel(statsG_temp)
%
for i = 1:size(statsG_temp,1)
% for i = 1:10
    disp(i)
    PixelList = statsG_temp(i).PixelList;
    PixelValues = statsG_temp(i).PixelValues;
    for j = 1:size(PixelList,1)
        y = PixelList(j,1);
        x = PixelList(j,2);
        z = PixelList(j,3);
        new_G(x,y,z) = 1;
    end
end
%%
%  imagesc(new_G(:,:,1))
% imwrite(new_G(:,:,1),[outpath '1.tif']);
%
new_G = logical(new_G);

F = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
F = logical(F);

size2 = Dis_thick;
parfor i = 1:num_images
    Dg = bwdistsc(new_G(:,:,i));
    curr = Dg<=(size2);
    curr = imfill(curr,4,'holes');
    F(:,:,i) = logical(uint8(curr) + uint8(BD(:,:,i)));
end
%
parfor i = 1:num_images
    BD_in = imfill(BD(:,:,i),'holes');
    all = ~BD_in;
    A = all | F(:,:,i);
    Area(i) = numel(find(~A));
end

%
for i = 1:num_images_2
    imwrite(double(F(:,:,i)),[outpath 'F_' sprintf('%03d',i) '.tif']);
    imwrite(logical(BD(:,:,i)),[exp_folder 'Result/0_BD/' 'BD_' sprintf('%03d',i) '.tif']);
end
save([exp_folder 'Result/0_BD/' 'Area.mat'],'Area');