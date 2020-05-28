clear all
%%
exp_folder = 'Z:\Jackie\neuroD6_retina_round2_section10_3\analysis\';
path = [exp_folder  'elastic_align/'];

stormpath = [path 'storm_merged/'];
stormpathds = [path 'storm_merged_ds/'];

stormfiles = [dir([stormpath '*.tif']) dir([stormpath '*.png'])];
stormfilesds = [dir([stormpathds '*.tif']) dir([stormpathds '*.png'])];
num_images = numel(stormfiles);
sinfo = imfinfo([stormpath stormfiles(1,1).name]);

outpath = 'Z:\Jackie\neuroD6_retina_round2_section10_3\9_4_2019_Processing\';

downsample = 10;
%set voxel size in nm
voxels=[(15.8), (15.8), 70];

BYs = zeros(ceil((sinfo(1,1).Height)/1),ceil((sinfo(1,1).Width)/1),num_images,'uint8');
disp('loading YFP')
%parpool(32)
parfor k = 1:num_images
    A = imread([stormpath stormfiles(k,1).name]);
    BYs(:,:,k) = A(:,:,3);
end
%%
figure;
imshow(max(BYs,[],3))
imwrite(max(BYs,[],3),[outpath 'yfp_storm_xy_all_raw_data.tif'])
imwrite(squeeze(max(BYs,[],1)),[outpath 'yfp_storm_xz_all_raw_data.tif'])
%%
gauss  = 5;
gausspix = (gauss);
sigmaZ = gausspix* voxels(1)/voxels(3);
sizZ= gausspix*2*sigmaZ* voxels(1)/voxels(3);
xZ=-ceil(sizZ/2):ceil(sizZ/2);
H1 = exp(-(xZ.^2/(2*sigmaZ^2)));
H1 = H1/sum(H1(:));
Hz=reshape(H1,[1 1 length(H1)]);
%create blurred gephyrin for region
%
by2 = zeros(size(BYs),'uint8');
parfor k=1:num_images
    disp(k)
    by2(:,:,k) = (imfilter(BYs(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
end
% tic

by2 = imfilter(by2,Hz,'same','replicate');
%%
figure;
imshow(max(by2,[],3))
%%
clear hy1
parfor j=1:num_images
    disp(j)
     A1 = by2(:,:,j); 
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
%%
figure;
hist(hy1dist,255)
%%
levelYs = multithresh(hy1dist,2);
leveluse = levelYs(2)/256
% leveluse = graythreshlist(uint8(hy1dist))

%
%parpool(32)
CYs = false(size(BYs));
disp('thresholding YFP')
parfor k = 1:num_images
    CYs(:,:,k) = im2bw(by2(:,:,k), leveluse);
end
%%
figure;
imshow(max(CYs,[],3))
%%
disp('conncomp YFP')
CCYs = bwconncomp(CYs,26);
disp('stats YFP')
%
statsYs1 = regionprops(CCYs,CYs,'Centroid', 'Area','PixelList','PixelIdxList','PixelValues');
%%
figure;
hist(log(1+[statsYs1.Area]),80)
%%
statsY2 = statsYs1(log(1+[statsYs1.Area])>10.5);

curr_stat = statsY2;
BY2 = zeros(size(BYs),'uint8');
disp('loading data')
for i = 1:numel(curr_stat)
    disp(int2str(i))
    BY2(curr_stat(i).PixelIdxList)=curr_stat(i).PixelValues*256;  
end

figure;
imshow((squeeze(max(BY2,[],3))))
%%
% statsY2_no = statsYs1(log(1+[statsYs1.Area])<=6);
% 
% curr_stat = statsY2_no;
% BY2_no = zeros(size(BYs),'uint8');
% disp('loading data')
% for i = 1:numel(curr_stat)
%     disp(int2str(i))
%     BY2_no(curr_stat(i).PixelIdxList)=curr_stat(i).PixelValues*256;  
% end
% 
% figure;
% imshow((squeeze(max(BY2_no,[],3))))
%%
%parpool(32)
%
clear temp
temp = statsYs1(log(1+[statsYs1.Area])<=10.5);
temp = temp(log(1+[temp.Area])>6);
clear B1
[B1(:,2),B1(:,1),B1(:,3)] = ind2sub(size(BY2),find(BY2));
%
B1_2 = B1(1:100:end,:);
%%
clear mindist 
for i=1:numel(temp)
    disp(i)
   mindisttemp = zeros(floor(numel([temp(i).PixelValues])/100),1);
    pixuse = zeros(floor(numel([temp(i).PixelValues])/100),3);
    for j=1:numel(mindisttemp)
        pixuse(j,:) = temp(i).PixelList(j*100,:);
    end
    parfor j=1:numel(mindisttemp)
    mindisttemp(j) = min(pdist2(pixuse(j,:),B1_2))
    end
   mindist(i) = min(mindisttemp);
end
%
figure;
hist(mindist,100)

%%
cut = 100;
statsY3 = cat(1,statsY2, temp(mindist<cut));
statsY4 = temp(mindist>=cut);
%
curr_stat = statsY3;
BY3 = zeros(size(BY2),'uint8');
disp('loading data')
for i = 1:numel(curr_stat)
    disp(int2str(i))
    BY3(curr_stat(i).PixelIdxList)=curr_stat(i).PixelValues*256;  
end

figure;
imshow(((max(BY3,[],3))))
%%
if exist([outpath 'allyfp/'])~=7
    mkdir([outpath 'allyfp/']);
    mkdir([outpath 'allyfp_small/']);
    mkdir([outpath 'nearconn_yfp/']);
    mkdir([outpath 'nearconn_yfp_small/']);
end
%
parfor k = 1:num_images
    disp(k)
    C3 = BY3(:,:,k);
    imwrite(C3, [outpath 'allyfp/' sprintf('%03d',k) '.tif']);
    C3_small = imresize(C3, 0.1);
    imwrite(C3_small, [outpath 'allyfp_small/' sprintf('%03d',k) '.tif']);
    
        C4 = BYs(:,:,k).*(BY3(:,:,k)/256);
    imwrite(C4, [outpath 'nearconn_yfp/' sprintf('%03d',k) '.tif']);
    C4_small = imresize(C4, 0.1);
    imwrite(C4_small, [outpath 'nearconn_yfp_small/' sprintf('%03d',k) '.tif']);
end