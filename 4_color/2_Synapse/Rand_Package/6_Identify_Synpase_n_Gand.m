%%
clear;clc;
outpath = 'Z:\Chenghang\chenghaz_014_XB2_P2_Control_B\analysis\Result\';
outpath2 = 'Z:\Chenghang\chenghaz_014_XB2_P2_Control_B\analysis\Result\Rand\';
load([outpath 'statsG2sw10.mat']);
statsGwater = statsGa2s;
%%
expfolder = 'Z:\Chenghang\chenghaz_014_XB2_P2_Control_B\analysis\elastic_align\storm_merged\';
files = [dir([expfolder '*.tif']) dir([expfolder '*.png'])];
infos = imfinfo([expfolder files(1,1).name]);
num_images = numel(files);
furtherds = 1;

new_G = zeros(ceil((infos(1,1).Height*furtherds)),ceil((infos(1,1).Width*furtherds)),num_images,'uint8');
%
frame = size(new_G);
statsGwater_random = statsGwater;
for i = 1:numel(statsGwater)
    disp(['Randomizing cluster # ' sprintf('%03d',i)]);
    PixelList = statsGwater(i).PixelList;
    Centroid = statsGwater(i).WeightedCentroid;
    mini = min(PixelList);
    maxi = max(PixelList);
    ran_x = randi(frame(2) - maxi(1) + mini(1));
    ran_y = randi(frame(1) - maxi(2) + mini(2));
    ran_z = randi(frame(3) - maxi(3) + mini(3));
    move_x = ran_x - mini(1);
    move_y = ran_y - mini(2);
    move_z = ran_z - mini(3);
    
    for j = 1:size(PixelList,1)
        statsGwater_random(i).PixelList(j,1) = PixelList(j,1) + move_x;
        statsGwater_random(i).PixelList(j,2) = PixelList(j,2) + move_y;
        statsGwater_random(i).PixelList(j,3) = PixelList(j,3) + move_z;
        statsGwater_random(i).PixelIdxList(j) = sub2ind(frame,statsGwater_random(i).PixelList(j,2),statsGwater_random(i).PixelList(j,1),statsGwater_random(i).PixelList(j,3));
    end
    
    statsGwater_random(i).WeightedCentroid = [Centroid(1)+move_x , Centroid(2)+move_y , Centroid(3)+move_z];
end

%
statsGwater_n_random = statsGwater_random;
save([outpath2 'statsG_n_rand.mat'],'statsGwater_n_random');