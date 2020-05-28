%%
%Read the filtered conventional images. 
expfolder = 'Z:\Jackie\neuroD6_retina_round2_section10_3\analysis\elastic_align\storm_merged\';
outpath = 'Z:\Jackie\neuroD6_retina_round2_section10_3\9_4_2019_Processing\temp\';
files = [dir([expfolder '*.tif']) dir([expfolder '*.png'])];
infos = imfinfo([expfolder files(1,1).name]);
num_images = numel(files);
furtherds = 1;

new_G = zeros(ceil((infos(1,1).Height*furtherds)),ceil((infos(1,1).Width*furtherds)),num_images,'uint8');

for i = 1:size(statsGa2ns,1)
% for i = 1:10
    disp(i)
    PixelList = statsGa2ns(i).PixelList;
    PixelValues = statsGa2ns(i).PixelValues;
    for j = 1:size(PixelList,1)
        y = PixelList(j,1);
        x = PixelList(j,2);
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
for i = 1:9
    imwrite(new_G(:,:,i),[outpath sprintf('%03d',i) '.tif']);
end