clear all
%%
base_path = 'Z:\Jackie\neuroD6_retina_round2_section10_3\';
exp_folder = [base_path 'analysis\'];

path = [exp_folder  'elastic_align\'];
mergedpath = [path 'storm_merged\'];
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);
disp(num_images)
%set voxel size in nm
voxel=[15.8, 15.8, 70];
outpath = 'Z:\Jackie\neuroD6_retina_round2_section10_3\9_4_2019_Processing\'
%%
clear BG
disp('allocating arrays')
BG = zeros(info.Height, info.Width, num_images,'uint8');
disp('loading data')
parfor j=1:num_images
    disp(j)
    AA = imread([mergedpath mergedfiles(j,1).name]);
    BG(:,:,j)= AA(:,:,2);
end
%%
gauss  = 5;
gausspix = (gauss);
sigmaZ = gausspix* voxel(1)/voxel(3);
sizZ= gausspix*2*sigmaZ* voxel(1)/voxel(3);
xZ=-ceil(sizZ/2):ceil(sizZ/2);
H1 = exp(-(xZ.^2/(2*sigmaZ^2)));
H1 = H1/sum(H1(:));
Hz=reshape(H1,[1 1 length(H1)]);
%create blurred gephyrin for region
%
bg2 = zeros(size(BG),'uint8');
parfor k=1:num_images
    disp(k)
    bg2(:,:,k) = (imfilter(BG(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
end

bg2 = imfilter(bg2,Hz,'same','replicate');
%%
figure;
imshow(imadjust(squeeze(max(bg2,[],3))));
%%
clear hy1
parfor j=1:num_images
    disp(j)
     A1 = bg2(:,:,j); 
     A1a = A1(find(A1));
     [hy, hx] = hist(A1a,0:1:255);
     hy1(j,:) = hy;
end
%
hx2 = 0:255;
hy1dist = [];
mhy1 = mean(hy1)./10;
for i=1:numel(mhy1)
nums = ones(1,round(mhy1(i))).*hx2(i);       
hy1dist = cat(2,hy1dist,nums);
end
%
figure;
hist(hy1dist,255)

threshfactorg = double(multithresh(hy1dist,2))
t_use = threshfactorg(2)/256;
%%
disp('making CG')
CG = false(size(bg2));
parfor k=1:size(bg2,3)
    CG(:,:,k) = im2bw(bg2(:,:,k), t_use);
end
%
%clear bg2
disp('making CCG')
CCG = bwconncomp(CG,26); 
%clear CG
disp('making statsG')
%
statsG = regionprops(CCG,BG,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
statsGgauss = regionprops(CCG,bg2,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
%%
%apply watershedding to statsG here


%volume_g0 = zeros(numel(statsG),1);
%volume_gall =zeros(numel(statsG),1);
%shape_g = zeros(numel(statsG),1);
%
%clear volume_g0 shape_g volume_gall
%for i=1:numel(statsG)
   % x = find([statsGgauss(jj).PixelValues]>1.2*threshfactorg(1));
   % volume_g0(jj) = numel(x);
    %volume_gall(jj) = numel(statsG(jj).PixelValues);
   % shape_g(jj) =numel(find([statsG(jj).PixelValues(x)>0])/volume_g0(jj));

%Volume2_t2b(i) = numel(find(statsGgauss(i,1).PixelValues>0.05*256));
%Volume1_t2b0(i) = numel(find(statsG(i,1).PixelValues(statsGgauss(i,1).PixelValues>...
%    0.05*256)>0));
%end
%val1w = (Volume1_t2b0(:)+1)./(Volume2_t2b(:)+1);
%val2w = log((Volume2_t2b(:)+1)*0.0158*0.0158*0.07);
%numel(find(val1w<0.99))
%figure;
%hist(val1w,100)
%
%Xn=80; Yn=80;
%Xrange=[0.3 1]; Yrange=[-8 -2];
%Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
%X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ; 

%figure;
%H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
%close
%
%cutoffg =200;
%H2 = H;
%H2(H2>cutoffg)=cutoffg;
%
%figure;
%pcolor(X,Y,H2)

%%
%this is used to see how many clusters fall above a certain size
figure;
hist(log([statsG.Area]+1),80)
%%
%manually examine clusters to set threshold to remove only cell bodies
areacutoff = multithresh(log(1+[statsG.Area]),1);
%%
%SHow each cluster
% for jj=1:numel(statsGgauss)
%     if log([statsG(jj).Area]+1)>areacutoff
%             disp(jj)
%             disp(log([statsG(jj).Area]+1))
%             minpix = min([statsGgauss(jj).PixelList]);
%             maxpix = max([statsGgauss(jj).PixelList]);
%             figure;
%             imshow(max(BG(minpix(2):maxpix(2),...
%                 minpix(1):maxpix(1),minpix(3):maxpix(3)),[],3))
%     end
% end
%%
statsGgauss = statsGgauss(log([statsG.Area]+1)>areacutoff);
statsG = statsG(log([statsG.Area]+1)>areacutoff);
%%
sizeBG = size(BG);
clear CG CCG
clear statsG2a

statsG2a(numel(statsGgauss)).stats = [];
clear minpix maxpix curr2 curr2b I2 I3 L L2 curr3 x y z

parfor jj=1:numel(statsGgauss)
   disp(jj)
minpix = min(statsGgauss(jj).PixelList);
maxpix = max(statsGgauss(jj).PixelList);
curr2 = zeros((maxpix(2)-minpix(2)+1+4),(maxpix(1)-minpix(1)+1+4),maxpix(3)- ...
minpix(3)+1+2,'uint8');
curr2b = zeros((maxpix(2)-minpix(2)+1+4),(maxpix(1)-minpix(1)+1+4),maxpix(3)- ...
minpix(3)+1+2,'uint8');
for j=1: numel(statsGgauss(jj).PixelList(:,1))
   curr2(statsG(jj).PixelList(j,2)-minpix(2)+1+2,statsG(jj).PixelList(j,1)-minpix(1)+1+2, ...
       statsG(jj).PixelList(j,3)-minpix(3)+1+1)=statsG(jj).PixelValues(j,1); 
   curr2b(statsGgauss(jj).PixelList(j,2)-minpix(2)+1+2,statsGgauss(jj).PixelList(j,1)-minpix(1)+1+2, ...
       statsGgauss(jj).PixelList(j,3)-minpix(3)+1+1)=statsGgauss(jj).PixelValues(j,1); 
end
I2 = imcomplement(curr2b);
I3 = imhmin(I2,10);
L = watershed(I3);
%clear currstat sizes
for i=1:max(L(:))
L2 = L;
curr3 = curr2b;
curr3(L~=i) = 0;
L2(L~=i) = 0;

[y, x, z] = ind2sub(size(curr3),find(curr3));
statsG2a(jj).stats(i).PixelList = [x y z];
for kk=1:numel(statsG2a(jj).stats(i).PixelList(:,1))
statsG2a(jj).stats(i).PixelValues(kk,1) = ...
    curr2(statsG2a(jj).stats(i).PixelList(kk,2), ...
    statsG2a(jj).stats(i).PixelList(kk,1),...
    statsG2a(jj).stats(i).PixelList(kk,3));
end

statsG2a(jj).stats(i).PixelList(:,1)= ...
    statsG2a(jj).stats(i).PixelList(:,1) + minpix(1) - 1-2;
statsG2a(jj).stats(i).PixelList(:,2)= ...
    statsG2a(jj).stats(i).PixelList(:,2) + minpix(2) - 1-2;
statsG2a(jj).stats(i).PixelList(:,3)= ...
    statsG2a(jj).stats(i).PixelList(:,3) + minpix(3) - 1-1;
statsG2a(jj).stats(i).Area = numel(statsG2a(jj).stats(i).PixelValues);
statsG2a(jj).stats(i).TintsG = sum(statsG2a(jj).stats(i).PixelValues);
statsG2a(jj).stats(i).PixelIdxList = ...
    sub2ind(sizeBG,statsG2a(jj).stats(i).PixelList(:,2),...
    statsG2a(jj).stats(i).PixelList(:,1),...
    statsG2a(jj).stats(i).PixelList(:,3));
end
end
%
clear templist
templist = []; statsG2a2 = []; cnt=0;
for i=1:numel(statsG2a)
templist = cat(1,templist, statsG2a(i).stats');
cnt = cnt+1;
     if cnt==2000
statsG2a2 = cat(1,statsG2a2, templist);
templist = [];
cnt = 0;
     end
end
statsG2a2 = cat(1,statsG2a2, templist);
statsGwater = statsG2a2;
%%
for i=1:numel(statsGwater)
    disp(i)
    statsGwater(i,1).PixelValues = BG(statsGwater(i,1).PixelIdxList);
    statsGwater(i,1).PixelValues2 = bg2(statsGwater(i,1).PixelIdxList);
statsGwater(i,1).Volume1_0 = numel(find(statsGwater(i,1).PixelValues>0));
statsGwater(i,1).Volume2_0 = numel(find(statsGwater(i,1).PixelValues2>0));
statsGwater(i,1).Volume2_t2a = numel(find(statsGwater(i,1).PixelValues2>(threshfactorg(1))));
statsGwater(i,1).Volume2_t2b = numel(find(statsGwater(i,1).PixelValues2>(1.1*threshfactorg(1))));
statsGwater(i,1).Volume2_t2c = numel(find(statsGwater(i,1).PixelValues2>(1.2*threshfactorg(1))));
statsGwater(i,1).Volume2_t2d = numel(find(statsGwater(i,1).PixelValues2>(1.3*threshfactorg(1))));
statsGwater(i,1).Volume2_t2e = numel(find(statsGwater(i,1).PixelValues2>(1.4*threshfactorg(1))));
statsGwater(i,1).Volume2_t2f = numel(find(statsGwater(i,1).PixelValues2>(1.5*threshfactorg(1))));
statsGwater(i,1).Volume2_t2g = numel(find(statsGwater(i,1).PixelValues2>(1.7*threshfactorg(1))));
statsGwater(i,1).Volume2_t2h = numel(find(statsGwater(i,1).PixelValues2>(2.0*threshfactorg(1))));
statsGwater(i,1).Volume1_t2a0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2b0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.1*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2c0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.2*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2d0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.3*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2e0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.4*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2f0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.5*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2g0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    1.7*threshfactorg(1))>0));
statsGwater(i,1).Volume1_t2h0 = ...
    numel(find(statsGwater(i,1).PixelValues(statsGwater(i,1).PixelValues2>...
    2.0*threshfactorg(1))>0));

statsGwater(i,1).Area = numel(statsGwater(i,1).PixelValues);
statsGwater(i,1).TintsG = sum(statsGwater(i,1).PixelValues);

statsGwater(i).WeightedCentroid(1) = ...
        sum([statsGwater(i).PixelList(:,1)].*...
        double([statsGwater(i).PixelValues]))/...
        (sum([statsGwater(i).PixelValues]));
statsGwater(i).WeightedCentroid(2) = ...
        sum([statsGwater(i).PixelList(:,2)].*...
        double([statsGwater(i).PixelValues]))/...
        (sum([statsGwater(i).PixelValues]));
statsGwater(i).WeightedCentroid(3) = ...
        sum([statsGwater(i).PixelList(:,3)].*...
        double([statsGwater(i).PixelValues]))/...
        (sum([statsGwater(i).PixelValues]));
end

numel(find([statsGwater.TintsG]>0))
statsGwater = statsGwater([statsGwater.TintsG]>0);
%
sizeshape_mat = cat(2,[statsGwater.Volume1_0]',...
    [statsGwater.Volume2_0]',[statsGwater.Volume2_t2a]',[statsGwater.Volume2_t2b]',...
    [statsGwater.Volume2_t2c]',[statsGwater.Volume2_t2d]',[statsGwater.Volume2_t2e]',...
    [statsGwater.Volume2_t2f]',[statsGwater.Volume2_t2g]',[statsGwater.Volume2_t2h]',...
    [statsGwater.Volume1_t2a0]', [statsGwater.Volume1_t2b0]',...
    [statsGwater.Volume1_t2c0]',[statsGwater.Volume1_t2d0]',...
    [statsGwater.Volume1_t2e0]',[statsGwater.Volume1_t2f0]',...
    [statsGwater.Volume1_t2g0]', [statsGwater.Volume1_t2h0]',...
    [statsGwater.Area]', [statsGwater.TintsG]');

centGw = zeros(numel(statsGwater),3);
dsvoxel = 158/0.3;
%dsvoxel = 158;
voxel = [15.8 15.8 70];
%
for jj =1:numel(statsGwater)   
    centGw(jj,:) = statsGwater(jj,1).WeightedCentroid;
end
%%
save([outpath 'sizeshapematG_and_cent_water10_area_cutoff.mat'],'centGw','sizeshape_mat')
save([outpath 'statsGwater10_area_cutoff.mat'],'statsGwater','-v7.3')
%%
clear all