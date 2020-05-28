%%
%some functions is required: (which should be under D:\Reaearch\Matlab\Storm_image_process\lib\)
    %Read_images
    %Make_mask
    %Z_resize
    %Rotate_Manually
%This file does: 
    %Read the raw images.
    %Make a mask on the images based on the intensity variation. (probably need modified)
    %Resize the Z-axis so that it has the same resolution as X-Y plane. 
    %(following codes are specific for this job)
    %Rotate the image manually to get a desire orientation. 
    %Delete pixels outside the on/off layers and rotation the plane back to
    %the original coordinate. 
%%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(28)
end
%%
%Load images
path='Z:\Chenghang\chenghaz005_2018.8.30\test\';
outpath = 'Z:\Chenghang\chenghaz005_2018.8.30\analysis\MATLAB\';
furtherds=1;
[BP,BG,BY] = Read_images(path,furtherds);
%%
%Make mask. Change the threshold (in percentage) in the Make_mask function.
num_images = size(BG,3);
[BG1,BY1]= Make_mask(BG,BY,num_images,0.005,0.005);
BP1 = BP;

%%
%remove bad sections. Bad sections can be determined manually and be
%removed here. 
badsec = [];

BP2 = BP1; BG2 = BG1; BY2 = BY1;
BP2(:,:,badsec) =[]; BG2(:,:,badsec) =[]; BY2(:,:,badsec) =[]; 
%%
% Resize on the z-axis. The default scale is *70/dsvoxel. 0.3 is the
% overall downsampling ratio. 
dsvoxel = 153/0.3;

[BP2a,BG2a,BY2a] = Z_resize(BP2,BG2,BY2,dsvoxel);
%%
%Rotate the whole image with three Eular angles. 0 <= angles < 360
xang = 4;
yang = 0;
zang = -4;
[BP3,BG3,BY3] = Rotate_Manually(xang,yang,zang,BP2a,BG2a,BY2a);

%%
%Further inspect all squeezed images after the rotation. Save the images.
%Be caureful for the path. Check if it need changed before running the
%code. 

clear C
C(:,:,1) = imadjust(squeeze(max(BP3,[],2)));
C(:,:,2) = imadjust(squeeze(max(BG3,[],2)));
C(:,:,3) = imadjust(squeeze(max(BY3,[],2)));
imwrite(C,[path 'rotated_ds_yz.tif'])
figure;
imshow(C)
%
clear C
C(:,:,1) = imadjust(squeeze(max(BP3,[],1)));
C(:,:,2) = imadjust(squeeze(max(BG3,[],1)));
C(:,:,3) = imadjust(squeeze(max(BY3,[],1)));
imwrite(C,[path 'rotated_ds_xz.tif'])
figure;
imshow(C)
%
clear C
C(:,:,1) = imadjust(squeeze(max(BP3,[],3)));
C(:,:,2) = imadjust(squeeze(max(BG3,[],3)));
C(:,:,3) = imadjust(squeeze(max(BY3,[],3)));
imwrite(C,[path 'rotated_ds_xy.tif'])
figure;
imshow(C)
%%
%manually determine the y-bound to get rid of some blank in images. 
ybounds_all = [6 86];%*_all is not currently used. 
ybounds_ipl = [16 76];
xbounds = [1 size(BY3,2)];
zbounds = [1 size(BY3,3)];
regionbounds_all = [ybounds_all(1) xbounds(1) zbounds(1);
                ybounds_all(2) xbounds(2) zbounds(2)];
regionbounds_ipl = [ybounds_ipl(1) xbounds(1) zbounds(1);
                ybounds_ipl(2) xbounds(2) zbounds(2)];

%%
%!!!probably wrong codes! Used B*3 rather than B*2 here!!!
BP4 = BP2; BG4 = BG2; BY4 = BY2;
V(1,:) = [regionbounds_ipl(1,1),regionbounds_ipl(1,2),regionbounds_ipl(1,3)];
V(2,:) = [regionbounds_ipl(1,1),regionbounds_ipl(2,2),regionbounds_ipl(1,3)];
V(3,:) = [regionbounds_ipl(1,1),regionbounds_ipl(1,2),regionbounds_ipl(2,3)];
V(4,:) = [regionbounds_ipl(2,1),regionbounds_ipl(1,2),regionbounds_ipl(1,3)];
V(5,:) = [regionbounds_ipl(2,1),regionbounds_ipl(2,2),regionbounds_ipl(1,3)];
V(6,:) = [regionbounds_ipl(2,1),regionbounds_ipl(1,2),regionbounds_ipl(2,3)];

[rV(:,2),rV(:,1),rV(:,3)]=tforminv(tform2,V(:,1),V(:,2),V(:,3));
rV(:,3) = rV(:,3)*(dsvoxel/70);
%
planepx1 = rV(1:3,1);   planepy1 = rV(1:3,2);   planepz1 = rV(1:3,3);
planepx2 = rV(4:6,1);   planepy2 = rV(4:6,2);   planepz2 = rV(4:6,3);
[xData1, yData1, zData1] = prepareSurfaceData( planepx1, planepy1, planepz1 );
[xData2, yData2, zData2] = prepareSurfaceData( planepx2, planepy2, planepz2 );
ft = fittype( 'poly11' );
[fitresult1] = fit( [xData1, yData1], zData1, ft );
[fitresult2] = fit( [xData2, yData2], zData2, ft );
coeffvals1 = coeffvalues(fitresult1);
coeffvals2 = coeffvalues(fitresult2);
clear currlist
[currlist(:,2), currlist(:,1),currlist(:,3)] = ind2sub(size(BP2),find(BP2));

out1 = currlist(:,3)>coeffvals1(1) + coeffvals1(2)*currlist(:,1) + ...
    coeffvals1(3)*currlist(:,2);
currgcl = currlist(out1,:);
out2 = currlist(:,3)<coeffvals2(1) + coeffvals2(2)*currlist(:,1) + ...
    coeffvals2(3)*currlist(:,2);
currinl = currlist(out2,:);

BP4(sub2ind(size(BP2),currgcl(:,2),currgcl(:,1),currgcl(:,3)))=0;
BP4(sub2ind(size(BP2),currinl(:,2),currinl(:,1),currinl(:,3)))=0;
BG4(sub2ind(size(BP2),currgcl(:,2),currgcl(:,1),currgcl(:,3)))=0;
BG4(sub2ind(size(BP2),currinl(:,2),currinl(:,1),currinl(:,3)))=0;
BY4(sub2ind(size(BP2),currgcl(:,2),currgcl(:,1),currgcl(:,3)))=0;
BY4(sub2ind(size(BP2),currinl(:,2),currinl(:,1),currinl(:,3)))=0;

%%
%Save your work! 
save([outpath 'rotated_ds_data.mat'],'badsec','tform2','BY2','BG2','BP2','BY3','BG3','outbounds','BP3','rot_tot',...
    'BY4','BG4','BP4','furtherds','dsvoxel');