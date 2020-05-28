clear all
%%
base_path = 'Z:\Jackie\neuroD6_retina_round2_section10_3\';
%base_path = 'Y:/backup/ysigal/';
exp_folder = [base_path 'analysis/'];
path = [exp_folder  'elastic_align/'];
stormpathds = [path 'storm_merged_ds/'];
convpathds = [path 'conv_merged_ds/'];
wgapathds = [path 'for_align_ds/'];
stormfilesds = [dir([stormpathds '*.tif']) dir([stormpathds '*.png'])];
convfilesds = [dir([convpathds '*.tif']) dir([convpathds '*.png'])];
wgafilesds = [dir([wgapathds '*.tif']) dir([wgapathds '*.png'])];
outpath = 'Z:\Jackie\neuroD6_retina_round2_section10_3\9_4_2019_Processing\';

infoc = imfinfo([convpathds convfilesds(1,1).name]);
infos = imfinfo([stormpathds stormfilesds(1,1).name]);
num_images = numel(stormfilesds);
%%
%load in and further downsample conventional images
furtherds = 0.3; dsvoxel = (158/furtherds);
BG = zeros(ceil((infoc(1,1).Height*furtherds)),ceil((infoc(1,1).Width*furtherds)),num_images,'uint8');
BP = zeros(ceil((infoc(1,1).Height*furtherds)),ceil((infoc(1,1).Width*furtherds)),num_images,'uint8');
BY = zeros(ceil((infoc(1,1).Height*furtherds)),ceil((infoc(1,1).Width*furtherds)),num_images,'uint8');
for k = 1:num_images 
    A = imread([convpathds convfilesds(k,1).name]);
    BP(:,:,k) = imresize(A(:,:,1),furtherds);
    BG(:,:,k) = imresize(A(:,:,2),furtherds);
    BY(:,:,k) = imresize(A(:,:,3),furtherds);
end
%%
%within image, determine if there is a local area of lower variance and
%make mask
BP1 = BP; BG1 = BG; BY1 = BY;
% for i=1:num_images
%     t1 = 60; t2 = 40;
% im1 = BP1(:,:,i);    im2 = BG1(:,:,i);
% filtim1 = (stdfilt(im1,ones(15,15)));  
% filtim2 = (stdfilt(im2,ones(15,15)));
% stdim(i) = numel(find(filtim1>t1))+ numel(find(filtim2>t2));
% if stdim(i)>50
%     filtuse = cat(1,(find(filtim2>t2)),(find(filtim1>t1)));
%     im1(filtuse)=0; im2(filtuse)=0;
%     BP1(:,:,i) = im1; BG1(:,:,i) = im2;
% figure;  filtusemat = zeros(size(im1),'uint8'); filtusemat(filtuse)=256;
% subplot(2,4,1); imshow(im2); subplot(2,4,2); imshow(mat2gray(filtim2))
% subplot(2,4,3); imshow((filtusemat));
% subplot(2,4,5);imshow(im1); subplot(2,4,6); imshow(mat2gray(filtim1))
% subplot(2,4,7); imshow((filtusemat));
% end
% end
%%
%remove bad sections
badsec = [66];
BP2 = BP1; BG2 = BG1; BY2 = BY1;
BP2(:,:,badsec) =[]; BG2(:,:,badsec) =[]; BY2(:,:,badsec) =[]; 
%%
Aa = permute(BP2,[1 3 2]);
Ab = permute(BG2,[1 3 2]);
Ac = permute(BY2,[1 3 2]);

for k = 1:ceil(size(Aa,3))
    A2a(:,:,k) = imresize(Aa(:,:,k), [ceil(size(Aa,1)) ...
        ceil(size(Aa,2)*furtherds*70/158)]);
    A2b(:,:,k) = imresize(Ab(:,:,k), [ceil(size(Ab,1)) ...
        ceil(size(Ab,2)*furtherds*70/158)]);
    A2c(:,:,k) = imresize(Ac(:,:,k), [ceil(size(Ac,1)) ...
        ceil(size(Ac,2)*furtherds*70/158)]);
end
BP2a = permute(A2a,[1 3 2]);
BG2a = permute(A2b,[1 3 2]);
BY2a = permute(A2c,[1 3 2]);
%%
figure;
imshow(max(BY2a,[],3))
%%
for    xang = 0;
    yang = 0;
    zang = -27;
xang = xang*pi/180;
yang = yang*pi/180;
zang = zang*pi/180;
clear Ypix
pcadim =2;
by3  =BY2a;
%[Ypix(:,1) Ypix(:,2) Ypix(:,3)] = ind2sub(size(by3),find(by3));

%[V2, D] = eig(cov(Ypix));
%D = diag(D);
%V = V2(:,:);

%v1 = [V(1,pcadim),V(2,pcadim),V(3,pcadim)];

%v2 = [0, 0, -1];
%
%r = vrrotvec(v1,v2);
%m = vrrotvec2mat(r)';
%
m= [cos(xang) 0 sin(xang) 0;
    0 1 0 0;
    -sin(xang) 0 cos(xang) 0
    0 0 0 1];
%
m1 = [1 0 0 0;
    0 cos(yang) -sin(yang) 0;
    0 sin(yang) cos(yang) 0;
    0 0 0 1];
m2 = [cos(zang) -sin(zang) 0 0;
    sin(zang) cos(zang) 0 0;
    0 0 1 0;
    0 0 0 1];
    %rot_tot = m1;
%tform = maketform('affine', rot_tot);

rot_tot = m1*m*m2;%transxy*m1;
tform = maketform('affine', rot_tot);
R = makeresampler('cubic','fill');
        
inbounds = [1 1 1;
            size(by3,1) size(by3,2) size(by3,3)];
outbounds = findbounds(tform,inbounds);
outysize = ceil(outbounds(2,1)-outbounds(1,1));
outxsize = ceil(outbounds(2,2)-outbounds(1,2));
outzsize = ceil(outbounds(2,3)-outbounds(1,3));
transxy2 = [1 0 0 0;
           0 1 0 0;
           0 0 1 0;
           -floor(outbounds(1,1))+1 -floor(outbounds(1,2))+1 -floor(outbounds(1,3))+1 1];
rot_tot2 = m1*m*m2*transxy2;
tform2 = maketform('affine', rot_tot2); 

BP3 = tformarray(BP2a,tform2,R,[1 2 3],[1 2 3],[outysize outxsize outzsize],[],0);
BG3 = tformarray(BG2a,tform2,R,[1 2 3],[1 2 3],[outysize outxsize outzsize],[],0);
BY3 = tformarray(BY2a,tform2,R,[1 2 3],[1 2 3],[outysize outxsize outzsize],[],0);

%%manually determine end of cell body

figure;
imshow(squeeze(max(BY3,[],2)))



clear C
C(:,:,1) = imadjust(squeeze(max(BP3,[],2)));
C(:,:,2) = imadjust(squeeze(max(BG3,[],2)));
C(:,:,3) = imadjust(squeeze(max(BY3,[],2)));
%imwrite(C,[path 'rotated_ds_yz.tif'])
figure;
imshow(C)
%
clear C
C(:,:,1) = imadjust(squeeze(max(BP3,[],1)));
C(:,:,2) = imadjust(squeeze(max(BG3,[],1)));
C(:,:,3) = imadjust(squeeze(max(BY3,[],1)));
%imwrite(C,[path 'rotated_ds_xz.tif'])
figure;
imshow(C)
%
clear C
C(:,:,1) = imadjust(squeeze(max(BP3,[],3)));
C(:,:,2) = imadjust(squeeze(max(BG3,[],3)));
C(:,:,3) = imadjust(squeeze(max(BY3,[],3)));
imwrite(C,[outpath 'rotated_ds_xy.tif'])
figure;
imshow(C)
end
%%
% ybounds_all = [6 86];
ybounds_ipl = [20 120];
xbounds = [1 size(BY3,2)];
zbounds = [1 size(BY3,3)];
% regionbounds_all = [ybounds_all(1) xbounds(1) zbounds(1);
%                 ybounds_all(2) xbounds(2) zbounds(2)];
regionbounds_ipl = [ybounds_ipl(1) xbounds(1) zbounds(1);
                ybounds_ipl(2) xbounds(2) zbounds(2)];
%%
BP4 = BP2; BG4 = BG2; BY4 = BY2;
V(1,:) = [regionbounds_ipl(1,1),regionbounds_ipl(1,2),regionbounds_ipl(1,3)];
V(2,:) = [regionbounds_ipl(1,1),regionbounds_ipl(2,2),regionbounds_ipl(1,3)];
V(3,:) = [regionbounds_ipl(1,1),regionbounds_ipl(1,2),regionbounds_ipl(2,3)];
V(4,:) = [regionbounds_ipl(2,1),regionbounds_ipl(1,2),regionbounds_ipl(1,3)];
V(5,:) = [regionbounds_ipl(2,1),regionbounds_ipl(2,2),regionbounds_ipl(1,3)];
V(6,:) = [regionbounds_ipl(2,1),regionbounds_ipl(1,2),regionbounds_ipl(2,3)];

[rV(:,2),rV(:,1),rV(:,3)]=tforminv(tform2,V(:,1),V(:,2),V(:,3));
rV(:,3) = rV(:,3)/furtherds*(158/70);
%
planepx1 = rV(1:3,1);   planepy1 = rV(1:3,2);   planepz1 = rV(1:3,3);
planepx2 = rV(4:6,1);   planepy2 = rV(4:6,2);   planepz2 = rV(4:6,3);
[xData1, yData1, zData1] = prepareSurfaceData( planepx1, planepy1, planepz1 );
[xData2, yData2, zData2] = prepareSurfaceData( planepx2, planepy2, planepz2 );
ft = fittype( 'poly11' );
[fitresult1, gof] = fit( [xData1, yData1], zData1, ft );
[fitresult2, gof] = fit( [xData2, yData2], zData2, ft );
coeffvals1 = coeffvalues(fitresult1);
coeffvals2 = coeffvalues(fitresult2);
clear currlist
[currlist(:,2), currlist(:,1) currlist(:,3)] = ind2sub(size(BP2),find(BP2));
%
out1 = currlist(:,3)>coeffvals1(1) + coeffvals1(2)*currlist(:,1) + ...
    coeffvals1(3)*currlist(:,2);
currgcl = currlist(out1,:);
out2 = currlist(:,3)<coeffvals2(1) + coeffvals2(2)*currlist(:,1) + ...
    coeffvals2(3)*currlist(:,2);
currinl = currlist(out2,:);
%
BP4(sub2ind(size(BP2),currgcl(:,2),currgcl(:,1),currgcl(:,3)))=0;
BP4(sub2ind(size(BP2),currinl(:,2),currinl(:,1),currinl(:,3)))=0;
BG4(sub2ind(size(BP2),currgcl(:,2),currgcl(:,1),currgcl(:,3)))=0;
BG4(sub2ind(size(BP2),currinl(:,2),currinl(:,1),currinl(:,3)))=0;
BY4(sub2ind(size(BP2),currgcl(:,2),currgcl(:,1),currgcl(:,3)))=0;
BY4(sub2ind(size(BP2),currinl(:,2),currinl(:,1),currinl(:,3)))=0;
%%            
save([outpath 'rotated_ds_data.mat'],'badsec','regionbounds_ipl',...
    'tform2','BY2','BG2','BP2','BY3','BG3','outbounds','BP3','rot_tot',...
    'BY4','BG4','BP4','furtherds','dsvoxel');