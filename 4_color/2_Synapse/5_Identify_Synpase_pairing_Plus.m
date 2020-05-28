clear;clc;
base_folder = 'Y:\Chenghang\04_4_Color\chenghaz_001_Sample_01_A\analysis\Result\';
outpath = base_folder;
%
load([outpath 'statsG2w10_edges_plus_2.mat'],'statsGwater')
tintsG_p140 = [statsGwater.tints_p140];
pairedg_idx = find(tintsG_p140);
numel(pairedg_idx)/numel(tintsG_p140)
statsGwater_sss = statsGwater(pairedg_idx);
save([base_folder 'G_paired_2.mat'],'statsGwater_sss');
%
load([outpath 'statsR2w10_edges_plus_2.mat'],'statsGwater')
tintsG_p140 = [statsGwater.tints_p140];
pairedg_idx = find(tintsG_p140);
numel(pairedg_idx)/numel(tintsG_p140)
statsRwater_sss = statsGwater(pairedg_idx);
save([base_folder 'R_paired_2.mat'],'statsRwater_sss');