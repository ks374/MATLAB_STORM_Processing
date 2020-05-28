%sections_chromalign
%%image_aligment%%
% define variables
%clear all
%arg1 = 'Y:/backup/ysigal/GLD3_bistratified_RGC/coverslips/cs16/';
%slicenum  = num2str(1);
local_exp =  arg1;
rel_conv_ints = '1111';
analysisfolder = cat(2, local_exp, 'analysis/');
ISanalysisfolder = cat(2, analysisfolder, 'individual_sections/');
acq_path = cat(2, local_exp, 'acquisition/')
slices = (numel(dir(fullfile(ISanalysisfolder, '0*')))-1);


conv_images_per_slice = numel(dir(fullfile([ISanalysisfolder ...
                        sprintf('%04d',str2num(slicenum)) '/rawimages/for_matlab/'],['*647IRconv_' sprintf('%03d',str2num(slicenum)) '_*'])))
    %
                    %
% load transformation matrix
reg_pth = [local_exp 'analysis/bead_fit/'];

load([reg_pth 'Visbeadworkspace_region2_10x']);

load([reg_pth 'IRbeadworkspace_region2_10x']);

%
for slice = str2num(slicenum);
for tile = 0:(conv_images_per_slice-1);
for pass = 0:1;
filename.storm488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/488storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filename.storm561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/561storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filename.storm647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/647storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filename.storm750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/750storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filename.convVis488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/488Visconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filename.convVis561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/561Visconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filename.convVis647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/647Visconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filename.convIR647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/647IRconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filename.convIR750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/rawimages/for_matlab/750IRconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filenameout.storm488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/488storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filenameout.storm647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/647storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filenameout.storm561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/561storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filenameout.storm750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/750storm_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_',sprintf('%01d',pass),'.tif');
filenameout.convVis488 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/488Visconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filenameout.convVis561 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/561Visconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filenameout.convVis647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/647Visconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filenameout.convIR647 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/647IRconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filenameout.convIR750 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/750IRconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_0.tif');
filenameout.convIR7502 = strcat(ISanalysisfolder,sprintf('%04d',slice),'/aligned/750IRconv_',sprintf('%03d',slice),'_',sprintf('%02d',tile),'_3.tif');

if exist(filenameout.convIR7502)==0
% load data
im.conv488 = im2double(imread(filename.convVis488));
im.conv488 = imresize(im.conv488,10)./ str2double(rel_conv_ints(4));
im.conv488adj = imadjust(im.conv488,stretchlim(im.conv488,[0 1]),[0 1]);

im.conv561 = im2double(imread(filename.convVis561));
im.conv561 = imresize(im.conv561,10)./str2double(rel_conv_ints(3));
im.conv561adj = imadjust(im.conv561,stretchlim(im.conv561,[0 1]),[0 1]);

im.conv647 = im2double(imread(filename.convVis647));
im.conv647 = imresize(im.conv647,10)./str2double(rel_conv_ints(2));
im.conv647adj = imadjust(im.conv647,stretchlim(im.conv647,[0 1]),[0 1]);

im.convIR647 = im2double(imread(filename.convIR647));
im.convIR647 = imresize(im.convIR647,10)./str2double(rel_conv_ints(2));
im.convIR647adj = imadjust(im.convIR647,stretchlim(im.convIR647,[0 1]),[0 1]);

im.conv750 = im2double(imread(filename.convIR750));
im.conv750 = imresize(im.conv750,10)./str2double(rel_conv_ints(1));
im.conv750adj = imadjust(im.conv750,stretchlim(im.conv750,[0 1]),[0 1]);
disp('images loaded')
if exist(filename.storm488)==2

im.storm488 = im2double(imread(filename.storm488));
else
im.storm488 = im.conv488adj;
end

%im.storm561 = im2double(imread(filename.storm561));
if exist(filename.storm647)==2
im.storm647 = im2double(imread(filename.storm647));
else
im.storm647 = im.conv647adj;
end

if exist(filename.storm750)==2
im.storm750 = im2double(imread(filename.storm750));
else
im.storm750 = im.conv750adj;
end

%correct for storm to conv drift
[output.storm488] = dftregistration(fft2(im.conv488adj),fft2(im.storm488),100);

%[output.storm561] = dftregistration(fft2(im.conv561adj),fft2(im.storm561),100);

[output.storm647] = dftregistration(fft2(im.convIR647adj),fft2(im.storm647),100);

[output.storm750] = dftregistration(fft2(im.conv750adj),fft2(im.storm750),100);

xform647 = [ 1  0  0
          0  1  0
         (output.storm647(4)) (output.storm647(3))  1 ];
tform_translate647 = maketform('affine',xform647);
xdata = [1 2560];
ydata = [1 2560];
imreg.storm647 = imtransform(im.storm647, tform_translate647, 'XData',xdata,'YData',ydata);

%xform561 = [ 1  0  0
%          0  1  0
%         (output.storm561(4)) (output.storm561(3))  1 ];
%tform_translate561 = maketform('affine',xform561);
%imreg.storm561 = imtransform(im.storm561, tform_translate561, 'XData',xdata,'YData',ydata);

xform488 = [  1  0  0
          0  1  0
         (output.storm488(4)) (output.storm488(3))  1 ];
tform_translate488 = maketform('affine',xform488);
[imreg.storm488] = imtransform(im.storm488, tform_translate488, 'XData',xdata,'YData',ydata);

xform750 = [  1  0  0
          0  1  0
        (output.storm750(4)) (output.storm750(3))  1 ];
tform_translate750 = maketform('affine',xform750);
[imreg.storm750] = imtransform(im.storm750, tform_translate750, 'XData',xdata,'YData',ydata);
disp('storm images aligned')
%pad images for warping
imtrans.conv488 = padarray(im.conv488, [2560 2560],'pre');
imtrans.conv561 = padarray(im.conv561, [0 2560],'pre');
imtrans.conv561 = padarray(imtrans.conv561, [2560 0],'post');
imtrans.conv750 = padarray(im.conv750, [0 2560],'post');
imtrans.conv750 = padarray(imtrans.conv750, [2560 0],'pre');

%if exist(filename.storm488)==2

imtrans.storm488 = padarray(imreg.storm488, [2560 2560],'pre');
%imtrans.storm561 = padarray(imreg.storm561, [0 2560],'pre');
%imtrans.storm561 = padarray(imtrans.storm561, [2560 0],'post');
imtrans.storm750 = padarray(imreg.storm750, [0 2560],'post');
imtrans.storm750 = padarray(imtrans.storm750, [2560 0],'pre');
%imshow(im.storm561)
%end

[im.conv488_warp] = imtransform(imtrans.conv488, tform_488_2_647,'XData',[1 2560],'YData',[1 2560]);
[im.conv561_warp] = imtransform(imtrans.conv561, tform_561_2_647,'XData',[1 2560],'YData',[1 2560]);
[im.conv750_warp] = imtransform(imtrans.conv750, tform_750_2_647,'XData',[1 2560],'YData',[1 2560]);

%if exist(filename.storm488)==2

[im.storm488_warp] = imtransform(imtrans.storm488, tform_488_2_647,'XData',[1 2560],'YData',[1 2560]);
%[im.storm561_warp] = imtransform(imtrans.storm561, tform_561_2_647,'XData',[1 2560],'YData',[1 2560]);
[im.storm750_warp] = imtransform(imtrans.storm750, tform_750_2_647,'XData',[1 2560],'YData',[1 2560]);
%end
%translate 647Visconv to 647IRconv
[output.conv647] = dftregistration(fft2(im.convIR647),fft2(im.conv647),100);
xform647 = [ 1  0  0
          0  1  0
         (output.conv647(4)) (output.conv647(3))  1 ];
tform_translate647 = maketform('affine',xform647);
xdata = [1 2560];
ydata = [1 2560];
im.conv647 = imtransform(im.conv647, tform_translate647, 'XData',xdata,'YData',ydata);
im.conv488_warp = imtransform(im.conv488_warp, tform_translate647, 'XData',xdata,'YData',ydata);
im.conv561_warp = imtransform(im.conv561_warp, tform_translate647, 'XData',xdata,'YData',ydata);

%if exist(filename.storm488)==2

im.storm488_warp = imtransform(im.storm488_warp, tform_translate647, 'XData',xdata,'YData',ydata);
%im.storm561_warp = imtransform(im.storm561_warp, tform_translate647, 'XData',xdata,'YData',ydata);
%end
disp('chromatic alignment done')
%save files
%im.rg_conv = cat(3, im.conv647, im.conv750_warp, im.conv488_warp);
%imwrite(im.rg_conv, filenameout.convrgb);
imwrite(im.conv488_warp, filenameout.convVis488);
imwrite(im.conv561_warp, filenameout.convVis561);
imwrite(im.conv647, filenameout.convVis647);
imwrite(im.convIR647, filenameout.convIR647);
imwrite(im.conv750_warp, filenameout.convIR750);
%if exist(filename.storm488)==2
imwrite(imreg.storm647, filenameout.storm647);
%imwrite(im.storm561_warp, filenameout.storm561);
imwrite(im.storm488_warp, filenameout.storm488);
imwrite(im.storm750_warp, filenameout.storm750);
%end


end
end
disp(['done writing slice', num2str(slice), num2str(tile)])
if exist('im') == 1
imwrite(imresize(im.conv750,0.05), filenameout.convIR7502);
end
end


end


