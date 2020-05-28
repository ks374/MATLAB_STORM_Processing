%%
%Read the filtered conventional images. 
expfolder = 'Y:\Chenghang\04_4_Color\chenghaz_001_Sample_01_A\analysis\elastic_align\storm_merged\';
outpath = 'Y:\Chenghang\04_4_Color\chenghaz_001_Sample_01_A\analysis\Result\5_V_Syn\';
files = [dir([expfolder '*.tif']) dir([expfolder '*.png'])];
infos = imfinfo([expfolder files(1,1).name]);
num_images = numel(files);
furtherds = 1;

new_G = zeros(ceil((infos(1,1).Height*furtherds)),ceil((infos(1,1).Width*furtherds)),num_images,'uint8');
%
statsG_temp = statsRwater_ssss;

for i = 1:size(statsG_temp,1)
% for i = 1:1000
    disp(i)
    PixelList = statsG_temp(i).PixelList;
    PixelValues = statsG_temp(i).PixelValues;
    for j = 1:size(PixelList,1)
        x = PixelList(j,2);
        y = PixelList(j,1);
        z = PixelList(j,3);
        new_G(x,y,z) = PixelValues(j);
    end
end
%     i = 1000;
%     PixelList = statsRwater(i).PixelList;
%     PixelValues = statsRwater(i).PixelValues;
%     for j = 1:size(PixelList,1)
%         y = PixelList(j,1);
%         x = PixelList(j,2);
%         z = PixelList(j,3);
%         new_G(x,y,z) = PixelValues(j);
%     end
%
% for i = 1:num_images
for i = 1:20
imwrite(new_G(:,:,i),[outpath 'R_' sprintf('%03d',i) '.tif']);
end