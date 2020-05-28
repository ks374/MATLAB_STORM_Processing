%%
%Read the filtered conventional images. 

outpath = 'Z:\Chenghang\0_For_Processing\chenghaz_011_Sample_03_P4_2\analysis\elastic_align\Result_Syn\';
statsG_temp = curr2;

num_images = numel(files);
new_G = zeros(ceil((infos(1,1).Height*furtherds)),ceil((infos(1,1).Width*furtherds)),num_images,'uint8');

for i = 1:size(statsG_temp,1)
% for i = 1
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
    Aa = permute(new_G,[1 3 2]);
    clear new_G
    for k = 1:ceil(size(Aa,3))
     A2a(:,:,k) = imresize(Aa(:,:,k), [ceil(size(Aa,1)) ...
         ceil(size(Aa,2)*70/15.5)]);
    end
    new_G = permute(A2a,[1 3 2]);
    
for i = 1:size(new_G,3)
%     i=1;
    imwrite(new_G(:,:,i),[outpath 'statsGwater_sss_' sprintf('%03d',i) '.tif']);
end