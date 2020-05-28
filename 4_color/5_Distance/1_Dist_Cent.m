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
%
load([outpath '5_V_Syn\R_paired_3.mat']);
statsRwater = statsRwater_ssss;
clear statsGwater
%
clear statsGa2s statsGwater
load([outpath '5_V_Syn/G_paired_3.mat']);
statsGwater = statsGwater_ssss;
%
WC_G = zeros(numel(statsGwater),3);
WC_R = zeros(numel(statsRwater),3);
for i = 1:numel(statsGwater)
    WC_G(i,:) = statsGwater(i).WeightedCentroid';
end
for i = 1:numel(statsRwater)
    WC_R(i,:) = statsRwater(i).WeightedCentroid';
end
%
WC_G(:,3) = WC_G(:,3)./70.*15.5;
WC_R(:,3) = WC_R(:,3)./70.*15.5;
%
DG = zeros(numel(statsGwater),1);
parfor i = 1:size(statsGwater,1)
    Z = squareform(pdist(cat(1,WC_G(i,:),WC_R)));
    Z1 = Z(:,1);
    Z1(1) = [];
    DG(i) = min(Z1);
end

DR = zeros(numel(statsRwater),1);
parfor i = 1:size(statsRwater,1)
    Z = squareform(pdist(cat(1,WC_R(i,:),WC_G)));
    Z1 = Z(:,1);
    Z1(1) = [];
    DR(i) = min(Z1);
end
%
if ~exist([outpath '6_Dist\'])
    mkdir([outpath '6_Dist\'])
end
save([outpath '6_Dist\Dist_Syn.mat'],'DG','DR');