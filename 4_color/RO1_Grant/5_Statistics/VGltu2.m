outpath = 'Y:\Chenghang\05_4_color_test\chenghaz_005_Sample_P8_C_A\analysis\elastic_align\Result\VGlut2\';
filepath = 'Y:\Chenghang\05_4_color_test\chenghaz_005_Sample_P8_C_A\analysis\elastic_align\storm_merged\';

files = [dir([filepath '*.tif']) dir([filepath '*.png'])];
infos = imfinfo([filepath files(1,1).name]);
num_images = numel(files);
furtherds = 1;
new_G = zeros(ceil((infos(1,1).Height*furtherds)),ceil((infos(1,1).Width*furtherds)),num_images,'uint8');
%
load([outpath 'V_paired']);
%
for i = 1:size(statsVwater_ss,1)
    disp(i)
    PixelList = statsVwater_ss(i).PixelList;
    PixelValues = statsVwater_ss(i).PixelValues;
    for j = 1:size(PixelList,1)
        x = PixelList(j,2);
        y = PixelList(j,1);
        z = PixelList(j,3);
        new_G(x,y,z) = PixelValues(j);
    end
end
%%
% imshow(new_G(:,:,1))
voxel = [15.5,15.5,70];
gauss  = 10;
gausspix = (gauss);
sigmaZ = gausspix* voxel(1)/voxel(3);
sizZ= gausspix*2*sigmaZ* voxel(1)/voxel(3);
xZ=-ceil(sizZ/2):ceil(sizZ/2);
H1 = exp(-(xZ.^2/(2*sigmaZ^2)));
H1 = H1/sum(H1(:));
Hz=reshape(H1,[1 1 length(H1)]);
%create blurred gephyrin for region
%
F = zeros(size(new_G),'uint8');
parfor k=1:num_images
    disp(k)
    F(:,:,k) = (imfilter(new_G(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
end

F = imfilter(F,Hz,'same','replicate');
%
CG = logical(F);
disp('making CCG')
CCG = bwconncomp(CG,26); 
%clear CG
disp('making statsG')
%
statsG = regionprops(CCG,new_G,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
%

for i = 1:numel(statsG)
    clear empty_pixel_idx
    empty_pixel_idx = find(statsG(i).PixelValues ==0);
    statsG(i).Area = statsG(i).Area - numel(empty_pixel_idx);
    statsG(i).PixelList(empty_pixel_idx,:) = [];
    statsG(i).PixelIdxList(empty_pixel_idx) = [];
    statsG(i).PixelValues(empty_pixel_idx) = [];
    statsG(i).TintsG = sum(statsG(i).PixelValues);
end
%%
statsVwater = statsG;
save([outpath 'Filtered_statsVwater.mat'] , 'statsVwater');