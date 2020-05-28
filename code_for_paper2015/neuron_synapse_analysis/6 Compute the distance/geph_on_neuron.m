clear all
%%
base_path = '/n/contefs1/backup/ysigal/';
%base_path = 'Y:/backup/ysigal/';

exp_folder = [base_path 'GLD3_bistratified_RGC/'];
%exp_folder = [base_path '03_22_13_sample1E_GABA_gephyrin/'];
%exp_folder = [base_path 'REIT_3rd_proc/'];
%exp_folder = [base_path 'REIT_3rd_proc/'];

%path = [exp_folder  'raw_data/'];
%path = [exp_folder  'cropped/'];
path = [exp_folder  'elastic_align/'];
path1 = [exp_folder  'elastic_thresh/'];
path2 = [path1  'parseyfp/'];
%
mergedpath = [path1 'storm_merged/'];
%
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);
disp(num_images)
%set voxel size in nm
voxel=[15.8, 15.8, 70];

%%
%load([path 'statsP2sw10_edges.mat'])
%
%statsPwater = statsPa2s;
%%
load([path1 'statsG2w10_edges_plus.mat'])
load([path1 'statsP2sw10_lo_thresh.mat'])
load([path1 'statslistG2sw10.mat']);
%%
load([path1 'statslistP2sw10_lo_thresh.mat']);
load([path1 'pairing_index_ps_gs_withedges.mat'],'pairedg_idx')
load([path1 'wga_normalization_pixels_withedges.mat'],'CW2')
%%
statsG2s = statsGwater;%(rcentGa2s(:,2)<46000);
Bfiles = dir([path2 'B1_' '*.mat']);
B1_2=cell(1,numel(Bfiles));
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
parfor i = 1:numel(Bfiles)
    disp(i)
    B1_2{1,i} = load(cat(2,path2,Bfiles(i).name))
end
   %
   B1_1=[];

voxel = [15.8 15.8 70];
for i=1:numel(Bfiles)
B1_1 = cat(1,B1_1,[B1_2{i}.x]);
end
B1_1a = B1_1.*repmat(voxel,numel(B1_1(:,1)),1);
%
figure;
hist(B1_1(:,2),100)
%
B2files = dir([path2 'B_' '*.mat']);
B2_2=cell(1,numel(B2files));
for i = 1:numel(B2files)
    disp(i)
    B2_2{1,i} = load(cat(2,path2,B2files(i).name));
end
B2_1=[];

for i=1:numel(B2files)
    disp(i)
B2_1 = cat(1,B2_1,[B2_2{i}.x]);
end
%
yfppath = [path1 'allyfp/'];
yfpfiles = [dir([yfppath '*.png']) dir([yfppath '*.tif'])];
num_images = numel(yfpfiles);
info = imfinfo([yfppath yfpfiles(1,1).name]);
BYs = zeros(info.Height,info.Width,num_images,'uint8');
Bidx = sub2ind(size(BYs),B2_1(:,2),B2_1(:,1),B2_1(:,3));
BYs(Bidx)=1;
voxel = [15.8 15.8 70];

%%
    dfiles = dir([path2 'D_*.mat']);
clear hd1 hd2
hxn1 = -1000:15.8:1000;
hxn2 = 0:20:4000;
hd1 = zeros(numel(dfiles),numel(hxn1));
hd2 = hd1; hd3 = hd1; hd4 = hd1;
hd5 = zeros(numel(dfiles),numel(hxn2));
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
parfor i=1:numel(dfiles)
disp(i)
x = load([path2 dfiles(i,1).name]);
    
%hd1(i,:) = histc(x.x(x.x(:,2)>1330 & x.x(:,2)<=2300,4),hxn1); %bottom
%hd2(i,:) = histc(x.x(x.x(:,2)<=1330,4),hxn1); %top
hd3(i,:) = histc(x.x(x.x(:,2)>100,4),hxn1); %no cellbody
hd4(i,:) = histc(x.x(:,4),hxn1); %all
hd5(i,:) = histc(x.x(:,3),hxn2);
end

%
shd5 = sum(hd5);
%
figure;
bar(hxn2,shd5,'histc')
%
shd1 = sum(hd1);
shd2 = sum(hd2);
shd3 = sum(hd3);
shd4 = sum(hd4);

%
figure;
bar(hxn1,shd3,'histc')
hold on
[xData1, yData1] = prepareCurveData( hxn1, shd3 );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 1e-05;
% Fit model to data.
[fitresult1, gof] = fit( xData1, yData1, ft, opts );
yfit1 = feval(fitresult1,xData1);
plot(xData1,yfit1)
shd = yfit1';
shd(shd==0)=1;
hd3a1 = shd.*(0.0158*0.0158*0.07);

figure;
bar(hxn1,shd4,'histc')
hold on
[xData1, yData1] = prepareCurveData( hxn1, shd4 );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 1e-05;
% Fit model to data.
[fitresult1, gof] = fit( xData1, yData1, ft, opts );
yfit1 = feval(fitresult1,xData1);
plot(xData1,yfit1)
shd = yfit1';
shd(shd==0)=1;
hd4a1 = shd.*(0.0158*0.0158*0.07);
%
%%

%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
parfor i=1:numel(statsG2s(:,1))
    disp(i)
    dist_y(i) = min(pdist2(statsG2s(i).WeightedCentroid.*voxel,B1_1a));
end
    idxnear_y = find(dist_y);
    %
    statsP2s = statsPa1s;%(rcentPa2s(:,2)<46000);
%
clear BP
disp('allocating arrays')

BP = zeros(info.Height, info.Width, num_images,'uint8');
disp('loading data')
for i = 1:numel(statsP2s)
    disp(int2str(i))
    BP(statsP2s(i).PixelIdxList)=statsP2s(i).PixelValues;  
end
%clear statsP2s statsPwater
    currstat = statsG2s(dist_y<1000);

    numel(currstat)
%
%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(currstat):-1:1
    disp(i)
currstat(i).PixelList_m140 = [];
currstat(i).RPixelValues_m140 = [];
currstat(i).YPixelValues_m140 = [];
currstat(i).wdistYR140 = [];
currstat(i).stdYR140 =  [];

currstat(i).YPixelValues = [];
currstat(i).wdistYG = [];
currstat(i).stdYG =  [];

currstat(i).PixelList_m70 = [];
currstat(i).RPixelValues_m70 = [];
currstat(i).YPixelValues_m70 = [];
currstat(i).wdistYR70 = [];
currstat(i).stdYR70 =  [];

%try preproccessing the BP varialbe in a for loop
   minpix = min(currstat(i).PixelList);  maxpix = max(currstat(i).PixelList);
   min1 = minpix(1)-100; min2 = minpix(2)-100; min3 = minpix(3)-15;
   max1 = maxpix(1)+100; max2 = maxpix(2)+100; max3 = maxpix(3)+15;
   if min1 < 1; min1=1; end 
   if min2 < 1; min2=1; end
   if min3 < 1; min3=1; end
   if max1 > info.Width; max1=info.Width; end
   if max2 > info.Height; max2=info.Height; end
   if max3 > num_images; max3=num_images; end      
   BPp(i).mat = BP(min2:max2,min1:max1,min3:max3); 
%try preproccessing the BY varialbe in a for loop    
   BYp(i).mat = BYs(min2:max2,min1:max1,min3:max3); 
end
%
%clear BYs BP
%

 for i=1:numel(BYp)
 test(i) = sum([BYp(i).mat(:)]);
 end
 %
 test2 = log(test+1);
 figure;
 hist(test2,30)
  for i=1:numel(BPp)
 test(i) = sum([BPp(i).mat(:)]);
 end
 %
 test2 = log(test+1);
 figure;
 hist(test2,30)

%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end   %
    parfor jj=1:numel(currstat)
        disp(jj)
        minpix = min(currstat(jj).PixelList); maxpix = max(currstat(jj).PixelList);
        min1 = minpix(1)-100; min2 = minpix(2)-100; min3 = minpix(3)-15;
        max1 = maxpix(1)+100; max2 = maxpix(2)+100; max3 = maxpix(3)+15;
        if min1 < 1; min1=1; end
        if min2 < 1; min2=1; end
        if min3 < 1; min3=1; end
        if max1 > info.Width; max1=info.Width; end
        if max2 > info.Height; max2=info.Height; end
        if max3 > num_images; max3=num_images; end        
    
    size1 = max1-min1 + 1; size2 = max2-min2 + 1; size3 = max3-min3 + 1;
            curr2 = false(size2, size1, size3);

        for j=1: numel(currstat(jj).PixelList(:,1))
            curr2(currstat(jj).PixelList(j,2)-min2+1, ...
               currstat(jj).PixelList(j,1)-min1+1, ...
               currstat(jj).PixelList(j,3)-min3+1)= 1;  
        end
        
          curr1a = BPp(jj).mat;
          curr1b = BYp(jj).mat;
   Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);

   size2= 140;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m140 = [x,y,z];
   
      size2= 70;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y, x, z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m70 = [x,y,z];
   
   RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
   for kk=1:numel(PixelList_m140(:,1))
      RPixelValues_m140(kk) = curr1a(PixelList_m140(kk,2),...
          PixelList_m140(kk,1), PixelList_m140(kk,3));
   end
   PixelList_m140(:,1)= PixelList_m140(:,1) + min1 - 1;
   PixelList_m140(:,2)= PixelList_m140(:,2) + min2- 1;
   PixelList_m140(:,3)= PixelList_m140(:,3) + min3 - 1;
   
    RPixelValues_m70 = zeros(numel(PixelList_m70(:,1)),1);
   for kk=1:numel(PixelList_m70(:,1))
      RPixelValues_m70(kk) = curr1a(PixelList_m70(kk,2),...
          PixelList_m70(kk,1), PixelList_m70(kk,3));
   end
   PixelList_m70(:,1)= PixelList_m70(:,1) + min1 - 1;
   PixelList_m70(:,2)= PixelList_m70(:,2) + min2- 1;
   PixelList_m70(:,3)= PixelList_m70(:,3) + min3 - 1;
        
    B1 = bwperim(curr1b,26);
    D = bwdistsc(B1,[voxel(1),voxel(2),voxel(3)]);
    D(curr1b>0)=-D(curr1b>0);
    
       YPixelValues = zeros(numel(currstat(jj).PixelList(:,1)),1);

   for kk=1:numel(currstat(jj).PixelList(:,1))
      YPixelValues(kk,1) =  D(currstat(jj).PixelList(kk,2)-min2+1, ...
                  currstat(jj).PixelList(kk,1)-min1+1, ...
                  currstat(jj).PixelList(kk,3)-min3+1);
   end
   currstat(jj).wdistYG = sum([YPixelValues].* ...
       double([currstat(jj).PixelValues]))/ ...
       (sum([currstat(jj).PixelValues]));
    currstat(jj).stdYG = std([YPixelValues(find([currstat(jj).PixelValues]>0))]);
       
       YPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);

    for kk=1:numel(PixelList_m140(:,1))
            YPixelValues_m140(kk,1) = ...
                D(PixelList_m140(kk,2)-min2+1, ...
                PixelList_m140(kk,1)-min1+1, ...
                PixelList_m140(kk,3)-min3+1);
    end
    
    YPixelValues_m70 = zeros(numel(PixelList_m70(:,1)),1);
     for kk=1:numel(PixelList_m70(:,1))
            YPixelValues_m70(kk,1) = ...
                D(PixelList_m70(kk,2)-min2+1, ...
                PixelList_m70(kk,1)-min1+1, ...
                PixelList_m70(kk,3)-min3+1);
     end
    
    currstat(jj).wdistYG = sum([YPixelValues].* ...
       double([currstat(jj).PixelValues]))/(sum([currstat(jj).PixelValues]));
   currstat(jj).stdYG =  std([YPixelValues(find([currstat(jj).PixelValues]>0))]);
     currstat(jj).YPixelValues = YPixelValues;

     
    currstat(jj).PixelList_m140 = PixelList_m140;
    currstat(jj).RPixelValues_m140 = RPixelValues_m140;
    currstat(jj).YPixelValues_m140 = YPixelValues_m140;
    currstat(jj).wdistYR140 = sum([YPixelValues_m140].* ...
       double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   currstat(jj).stdYR140 =  std([YPixelValues_m140(find([RPixelValues_m140]>0))]);
   
       currstat(jj).PixelList_m70 = PixelList_m70;
    currstat(jj).RPixelValues_m70 = RPixelValues_m70;
    currstat(jj).YPixelValues_m70 = YPixelValues_m70;
    currstat(jj).wdistYR70 = sum([YPixelValues_m70].* ...
       double([RPixelValues_m70]))/(sum([RPixelValues_m70]));
   currstat(jj).stdYR70 =  std([YPixelValues_m70(find([RPixelValues_m70]>0))]);
    end
    %
        statsGs_near = currstat;

     save([path 'statsG_nearw10.mat'],...
     'statsGs_near','-v7.3')
 %
save([path 'disty_for_statsG2w10_plus.mat'],...
     'dist_y','-v7.3')

%%
figure;
hist(dist_y,1000)
%%
paired_near = pairedg_idx(dist_y<1000);
statsGs_near_p = statsGs_near(paired_near);
statsGs_near_up = statsGs_near(~paired_near);
%%
hd2c1 = histc([statsGs_near_p.wdistYG],hxn1);
hd3c1 = hd2c1./hd3a1;
hd2c2 = histc([statsGs_near_p.wdistYR70],hxn1);
hd3c2 = hd2c2./hd3a1;
hd2c3 = histc([statsGs_near_up.wdistYG],hxn1);
hd3c3 = hd2c3./hd3a1;
figure;
plot(hxn1,(hd3c1),'g')
hold on
plot(hxn1,(hd3c2),'r')
hold on
plot(hxn1,(hd3c3),'b')
figure;
plot(hxn1,smooth(hd3c1,3),'g')
hold on
plot(hxn1,smooth(hd3c2,3),'r')
hold on
plot(hxn1,smooth(hd3c3,3),'b')
savefig([path 'density_distY_paired_g(g)_and_p70(r)and_unpaired_g(b).fig'])
%%

%%
lowbound = -140;
hibound = 70;
statsGon_p = statsGs_near_p([statsGs_near_p.wdistYG]>lowbound & ...
    [statsGs_near_p.wdistYG]<hibound);
statsGon_up = statsGs_near_up([statsGs_near_up.wdistYG]>lowbound & ...
    [statsGs_near_up.wdistYG]<hibound);
statsGnearonly = statsGs_near_p([statsGs_near_p.wdistYG]>hibound);
numel(statsGon_p)/(numel(statsGon_up)+ numel(statsGon_p))
numel(statsGon_p)

%%
hd2c1 = histc([statsGon_p.wdistYG],hxn1);
hd3c1 = hd2c1./hd3a1;
hd2c2 = histc([statsGon_p.wdistYR140],hxn1);
hd3c2 = hd2c2./hd3a1;
figure;
plot(hxn1,(hd3c1),'g')
hold on
plot(hxn1,(hd3c2),'r')
%savefig([path 'density_distY_on_paired_g(g)_and_p140(r).fig'])
%%
figure;
[hy hx] = hist([statsGon_p.wdistYR140]-[statsGon_p.wdistYG],-250:20:250);
hold on
[hy2 hx] = hist([statsGnearonly.wdistYR140]-[statsGnearonly.wdistYG],-250:20:250);
plot(hx,hy./sum(hy))
hold on
plot(hx,hy2./sum(hy2))
savefig([path 'delta_distY_on_paired_and_near.fig'])
%%

save([path 'statsGon_final.mat'],'statsGon_p','statsGon_up')
