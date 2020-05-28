% this script is used to generate a bead map and corresponding warping 
% transforms based on beads collected with the quadview and visible in
% 750and 647 channels (100nm 715/755 beads).  Transforms map 750 to the 647 channel and 
%  .bin files are expected in the folder "region1"
% and the xml file from dave used to collect the beads should also be named 
% 'region1' and is expected in the same folder. File names should contain 'IR'
% If transforms will be applied to images with pixel size other than 256x256 pixels (158nm/pixel) ,
% adjust the scale accordingly.
%user parameters
%clear all
%arg1 = 'Y:/ysigal/sample1B_glyR_geph/coverslips/sample1B_glyR_geph_cs10/';

local_exp =  arg1;
acq_path = [local_exp 'acquisition/'];
base_pth = [local_exp 'analysis/bins/'];
save_pth = [local_exp 'analysis/bead_fit/'];
scale = 10;

xdiv = 256; % width of left channel in pixels
ydiv = 256; % height of left channel in pixels
match_radius = 6*scale; % radius tolerance for matching left/right channels
init_750xshift = 0; % pixel shift between 647 and 750 channel
init_750yshift = 256*scale; % pixel shift between 647 and 750 channels
init_561xshift = 256*scale; % pixel shift between 647 and 561 channel
init_561yshift = 0; % pixel shift between 647 and 561 channels
init_488xshift = 256*scale; % pixel shift between 647 and 488 channel
init_488yshift = 256*scale; % pixel shift between 647 and 488 channels
nm_per_pix = 158/scale; % nanometers per pixel

ax_fs = 10; % axis fontsize for plots
la_fs = 10; % label fontsize for plots
ti_fs = 10; % title fontsize for plots

%%IRtransforms%%
% file i/o stuff
numIRbeadmov = numel(dir([acq_path 'IRbeads*' '*.dax']));
numVisbeadmov = numel(dir([acq_path 'Visbeads*' '*.dax']));
if numVisbeadmov>1
Vis_workspace_nm_r2 = [save_pth 'Visbeadworkspace_region2_' num2str(scale) 'x.mat'];
load(Vis_workspace_nm_r2);
tform_561start = tform_561_2_647
tform_488start = tform_488_2_647
end

if numIRbeadmov>1
IR_workspace_nm_r2 = [save_pth 'IRbeadworkspace_region2_' num2str(scale) 'x.mat'];    
load(IR_workspace_nm_r2);
tform_750start = tform_750_2_647


match_radius = 0.5*scale;


if exist([save_pth '750_2_647_total_offsets.tif'],'file')==0

track_file_nms = dir([base_pth 'IRbeads*' '_1_' '*mlist.bin'] );
len = length(track_file_nms);

for k=1:len
track_file_nm = [base_pth track_file_nms(k).name];
    mol_list_out_nm = [track_file_nm(1:end-4) '.mat'];

 % load molecule list from insight3
    %mol_list_out_nm = [mol_list_out_nm(1:end-4) '.mat'];
    M = bin2mat(track_file_nm,mol_list_out_nm,0); % thanks wenqin
    A(k).M = M;
    A(k).bin_nm = track_file_nm;
    
    % keep only those molecules which are long enough
    keep_ind = [];
    min_trail_lth = (length(A(k).M(1).x)-1); % minimum length of trail for localizations to be used

    for m = 1:length(M)
        lth = length(M(m).x)-1;
        if lth>=min_trail_lth
            keep_ind = [keep_ind m];
        end
    end
     A(k).M = M(keep_ind);
     
     
    A(k).M2 = (A(k).M); % average of bead position during trace

    for b = 1:length(A(k).M)
        A(k).M2(1,b).x = (A(k).M2(1,b).x*scale);
        A(k).M2(1,b).y = (A(k).M2(1,b).y*scale);
        
    end
    
    A(k).set647_inds = find([A(k).M2.xc]>0 & [A(k).M2.xc]<=xdiv & [A(k).M2.yc]<=ydiv);
    A(k).set647_pos.x = double([A(k).M2(A(k).set647_inds).x]);
    A(k).set647_pos.y = double([A(k).M2(A(k).set647_inds).y]);
    
    A(k).set750_inds = find([A(k).M2.xc]>0 & [A(k).M2.xc]<=xdiv & [A(k).M2.yc]>=ydiv);
    A(k).set750_pos.x = double([A(k).M2(A(k).set750_inds).x]);
    A(k).set750_pos.y = double([A(k).M2(A(k).set750_inds).y]);

    [A(k).matched750,A(k).unmatched750] = corr_mols(A(k).set750_pos,A(k).set647_pos,tform_750start,match_radius);

    
disp([num2str(2*length(A(k).matched750.set1_inds)) '/' ...
        num2str(2*length(A(k).matched750.set1_inds) ...
        +length(A(k).unmatched750.set1_inds)+ ...
        length(A(k).unmatched750.set2_inds)) ' 750 molecules matched'])
        disp(['bead image number' num2str(k)])
        
    A(k).S750_1 = A(k).M2(A(k).set750_inds(A(k).matched750.set1_inds)); % matched molecules, set1
    A(k).S750_2 = A(k).M2(A(k).set647_inds(A(k).matched750.set2_inds)); % matched molecules, set2
end

%%% group together into two big sets
comb_set750_1_pos.x = []; % set1 = left channel (pixel units)
comb_set750_1_pos.y = []; 
comb_set750_2_pos.x = []; % set2 = right channel (pixel units)
comb_set750_2_pos.y = [];


for k = 1:length(A)
    try
        comb_set750_1_pos.x = double([comb_set750_1_pos.x A(k).S750_1(:).x]);
        comb_set750_1_pos.y = double([comb_set750_1_pos.y A(k).S750_1(:).y]);
        comb_set750_2_pos.x = double([comb_set750_2_pos.x A(k).S750_2(:).x]);
        comb_set750_2_pos.y = double([comb_set750_2_pos.y A(k).S750_2(:).y]);
    catch
    end
end

% find mapping from set2 (right) onto set1 (left), save tform_right2left
set750_1_points = [comb_set750_1_pos.x' comb_set750_1_pos.y'];
set750_2_points = [comb_set750_2_pos.x' comb_set750_2_pos.y'];

[warped_set750_2_pos.x,warped_set750_2_pos.y] = tforminv(tform_750_2_647,comb_set750_2_pos.x,comb_set750_2_pos.y);

set750_2_warp_error_x = comb_set750_1_pos.x - warped_set750_2_pos.x;
set750_2_warp_error_y = comb_set750_1_pos.y - warped_set750_2_pos.y;
std_set750_2_warp_error_x = std(set750_2_warp_error_x);
std_set750_2_warp_error_y = std(set750_2_warp_error_y);

set750_2_orig_error_x = comb_set750_1_pos.x - comb_set750_2_pos.x - init_750xshift;
set750_2_orig_error_y = comb_set750_1_pos.y - comb_set750_2_pos.y - init_750yshift;
std_set750_2_orig_error_x = std(set750_2_orig_error_x);
std_set750_2_orig_error_y = std(set750_2_orig_error_y);


std_dist = nm_per_pix*std(sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2));
mean_dist = nm_per_pix*mean(sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2));


[cdf750.y, cdf750.x] = ecdf(nm_per_pix*sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2));
cdf90_750 = (cdf750.x(find(cdf750.y>0.9,1,'first')));
disp(['90% of 750 beads aligned to ' num2str(cdf90_750) 'nm ,using ' num2str(length(set750_1_points(:,1))) ' beads'])


% plot deviations before warping for right to left, stage to left registration
figure(1); clf
set(gcf,'position',[493 441 1112 482])
fac = 2;
hold on
plot(comb_set750_2_pos.x,comb_set750_2_pos.y,'g.')
quiver(comb_set750_2_pos.x,comb_set750_2_pos.y,fac*set750_2_orig_error_x,fac*set750_2_orig_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['750 to 647, no warping ' num2str(fac) 'x mag. residuals' ],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
saveas(gcf,[save_pth 'IRbead_errors_pre_warping.fig'],'fig')
saveas(gcf,[save_pth 'IRbead_errors_pre_warping.tif'],'tif')

% plot deviations after warping for right to left, stage to left registration
figure(2); clf
set(gcf,'position',[493 441 1212 482])
fac = 10;
hold on
plot(comb_set750_2_pos.x,comb_set750_2_pos.y,'g.')
quiver(comb_set750_2_pos.x,comb_set750_2_pos.y,fac*set750_2_warp_error_x,fac*set750_2_warp_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['750 to 647  ' num2str(fac) ...
    'x mag. residuals'],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
saveas(gcf,[save_pth 'IRbead_errors_post_warping.fig'],'fig')
saveas(gcf,[save_pth 'IRbead_errors_post_warping.tif'],'tif')

% plot x and y offsets for right to left warping
figure(3); clf
set(gcf,'position',[493 441 1112 482])
subplot(1,2,1);
hist(set750_2_warp_error_x*nm_per_pix,15)
set(gca,'fontsize',ax_fs)
xlabel('measured x offsets [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title(['750 to 647 STD_x=' num2str(std_set750_2_warp_error_x*nm_per_pix,'%0.1f') ' nm'],'fontsize',ti_fs)
set(gca,'position',[0.1544 0.2 0.3102 0.65])
subplot(1,2,2);
hist(set750_2_warp_error_y*nm_per_pix,15)
set(gca,'fontsize',ax_fs)
xlabel('measured y offsets [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title(['750 to 647 STD_y=' num2str(std_set750_2_warp_error_y*nm_per_pix,'%0.1f') ' nm'],'fontsize',ti_fs)
set(gca,'position',[0.5948 0.2 0.3102 0.65])
saveas(gcf,[save_pth '750_2_647_xy_offsets.fig'],'fig')
saveas(gcf,[save_pth '750_2_647_xy_offsets.tif'],'tif')

% plot total offsets for left to right warping
figure(4); clf
set(gcf,'position',[493 441 1112 482])
subplot(1,2,1)
hist(nm_per_pix*sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2),40)
set(gca,'fontsize',ax_fs)
xlabel('measured displacements [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title('750 to 647 measured displacements','fontsize',ti_fs)
subplot(1,2,2)
cdfplot(nm_per_pix*sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2))
set(gca,'fontsize',ax_fs)
xlabel('distance error [nm]','fontsize',la_fs)
ylabel('cumulative probability','fontsize',la_fs)
title('Cumulative distribution','fontsize',la_fs)
saveas(gcf,[save_pth '750_2_647_total_offsets.fig'],'fig')
saveas(gcf,[save_pth '750_2_647_total_offsets.tif'],'tif')
end
end

if numVisbeadmov>1
if exist([save_pth '488_2_647_total_offsets.tif'],'file')==0

track_file_nms = dir([base_pth 'Visbeads*' '_1_' '*mlist.bin'] );
len = length(track_file_nms);

for k=1:len
track_file_nm = [base_pth track_file_nms(k).name];
    mol_list_out_nm = [track_file_nm(1:end-4) '.mat'];

 % load molecule list from insight3
    M = bin2mat(track_file_nm,mol_list_out_nm,0); % thanks wenqin
    A(k).M = M;
    A(k).bin_nm = track_file_nm;
    
    % keep only those molecules which are long enough
    keep_ind = [];
    min_trail_lth = (length(A(k).M(1).x)-1); % minimum length of trail for localizations to be used

    for m = 1:length(M)
        lth = length(M(m).x)-1;
        if lth>=min_trail_lth
            keep_ind = [keep_ind m];
        end
    end
     A(k).M = M(keep_ind);
     
     
    A(k).M2 = (A(k).M); % average of bead position during trace

    for b = 1:length(A(k).M)
        A(k).M2(1,b).x = (A(k).M2(1,b).x*scale);
        A(k).M2(1,b).y = (A(k).M2(1,b).y*scale);
        
    end

    
    A(k).set647_inds = find([A(k).M2.xc]>0 & [A(k).M2.xc]<=xdiv & [A(k).M2.yc]<=ydiv);
    A(k).set647_pos.x = double([A(k).M2(A(k).set647_inds).x]);
    A(k).set647_pos.y = double([A(k).M2(A(k).set647_inds).y]);
    
    A(k).set561_inds = find([A(k).M2.xc]>0 & [A(k).M2.xc]>=xdiv & [A(k).M2.yc]<=ydiv);
    A(k).set561_pos.x = double([A(k).M2(A(k).set561_inds).x]);
    A(k).set561_pos.y = double([A(k).M2(A(k).set561_inds).y]);
    
    A(k).set488_inds = find([A(k).M2.xc]>0 & [A(k).M2.xc]>=xdiv & [A(k).M2.yc]>=ydiv);
    A(k).set488_pos.x = double([A(k).M2(A(k).set488_inds).x]);
    A(k).set488_pos.y = double([A(k).M2(A(k).set488_inds).y]);


 [A(k).matched561,A(k).unmatched561] = corr_mols(A(k).set561_pos,A(k).set647_pos,tform_561start,match_radius);
  [aa ab ac] = unique(A(k).matched561.set2_inds);
 A(k).matched561.set1_inds = A(k).matched561.set1_inds(ab);
 A(k).matched561.set2_inds = A(k).matched561.set2_inds(ab); 
 [A(k).matched488,A(k).unmatched488] = corr_mols(A(k).set488_pos,A(k).set647_pos,tform_488start,match_radius);
   [aa ab ac] = unique(A(k).matched488.set2_inds);
 A(k).matched488.set1_inds = A(k).matched488.set1_inds(ab);
 A(k).matched488.set2_inds = A(k).matched488.set2_inds(ab);
    
  A(k).S561_1 = A(k).M2(A(k).set561_inds(A(k).matched561.set1_inds)); % matched molecules, set1
  A(k).S561_2 = A(k).M2(A(k).set647_inds(A(k).matched561.set2_inds)); % matched molecules, set2
  
  A(k).S488_1 = A(k).M2(A(k).set488_inds(A(k).matched488.set1_inds)); % matched molecules, set1
  A(k).S488_2 = A(k).M2(A(k).set647_inds(A(k).matched488.set2_inds)); % matched molecules, set2
    
end


%%% group together into two big sets
comb_set561_1_pos.x = []; % set1 = left channel (pixel units)
comb_set561_1_pos.y = []; 
comb_set561_2_pos.x = []; % set2 = right channel (pixel units)
comb_set561_2_pos.y = [];

comb_set488_1_pos.x = []; % set1 = left channel (pixel units)
comb_set488_1_pos.y = []; 
comb_set488_2_pos.x = []; % set2 = right channel (pixel units)
comb_set488_2_pos.y = [];


for k = 1:length(A)
    try
        comb_set561_1_pos.x = double([comb_set561_1_pos.x A(k).S561_1(:).x]);
        comb_set561_1_pos.y = double([comb_set561_1_pos.y A(k).S561_1(:).y]);
        comb_set561_2_pos.x = double([comb_set561_2_pos.x A(k).S561_2(:).x]);
        comb_set561_2_pos.y = double([comb_set561_2_pos.y A(k).S561_2(:).y]);
     
        comb_set488_1_pos.x = double([comb_set488_1_pos.x A(k).S488_1(:).x]);
        comb_set488_1_pos.y = double([comb_set488_1_pos.y A(k).S488_1(:).y]);
        comb_set488_2_pos.x = double([comb_set488_2_pos.x A(k).S488_2(:).x]);
        comb_set488_2_pos.y = double([comb_set488_2_pos.y A(k).S488_2(:).y]);
    catch
    end
end

% find mapping from set2 (right) onto set1 (left), save tform_right2left
set561_1_points = [comb_set561_1_pos.x' comb_set561_1_pos.y'];
set561_2_points = [comb_set561_2_pos.x' comb_set561_2_pos.y'];

set488_1_points = [comb_set488_1_pos.x' comb_set488_1_pos.y'];
set488_2_points = [comb_set488_2_pos.x' comb_set488_2_pos.y'];



[warped_set561_2_pos.x,warped_set561_2_pos.y] = tforminv(tform_561_2_647,comb_set561_2_pos.x,comb_set561_2_pos.y);

set561_2_warp_error_x = comb_set561_1_pos.x - warped_set561_2_pos.x;
set561_2_warp_error_y = comb_set561_1_pos.y - warped_set561_2_pos.y;
std_set561_2_warp_error_x = std(set561_2_warp_error_x);
std_set561_2_warp_error_y = std(set561_2_warp_error_y);

set561_2_orig_error_x = comb_set561_1_pos.x - comb_set561_2_pos.x - init_561xshift;
set561_2_orig_error_y = comb_set561_1_pos.y - comb_set561_2_pos.y - init_561yshift;
std_set561_2_orig_error_x = std(set561_2_orig_error_x);
std_set561_2_orig_error_y = std(set561_2_orig_error_y);


[warped_set488_2_pos.x,warped_set488_2_pos.y] = tforminv(tform_488_2_647,comb_set488_2_pos.x,comb_set488_2_pos.y);

set488_2_warp_error_x = comb_set488_1_pos.x - warped_set488_2_pos.x;
set488_2_warp_error_y = comb_set488_1_pos.y - warped_set488_2_pos.y;
std_set488_2_warp_error_x = std(set488_2_warp_error_x);
std_set488_2_warp_error_y = std(set488_2_warp_error_y);

set488_2_orig_error_x = comb_set488_1_pos.x - comb_set488_2_pos.x - init_488xshift;
set488_2_orig_error_y = comb_set488_1_pos.y - comb_set488_2_pos.y - init_488yshift;
std_set488_2_orig_error_x = std(set488_2_orig_error_x);
std_set488_2_orig_error_y = std(set488_2_orig_error_y);

[cdf488.y, cdf488.x] = ecdf(nm_per_pix*sqrt(set488_2_warp_error_x.^2 + set488_2_warp_error_y.^2));
cdf90_488 = (cdf488.x(find(cdf488.y>0.9,1,'first')));
disp(['90% of 488 beads aligned to ' num2str(cdf90_488) 'nm ,using ' num2str(length(set488_1_points(:,1))) ' beads'])

[cdf561.y, cdf561.x] = ecdf(nm_per_pix*sqrt(set561_2_warp_error_x.^2 + set561_2_warp_error_y.^2));
cdf90_561 = (cdf561.x(find(cdf561.y>0.9,1,'first')));
disp(['90% of 561 beads aligned to ' num2str(cdf90_561) 'nm ,using ' num2str(length(set488_1_points(:,1))) ' beads'])


% plot deviations before warping for right to left, stage to left registration
figure(1); clf
set(gcf,'position',[493 441 1112 482])
fac = 2;
subplot(1,2,1)
hold on
plot(comb_set561_2_pos.x,comb_set561_2_pos.y,'g.')
quiver(comb_set561_2_pos.x,comb_set561_2_pos.y,fac*set561_2_orig_error_x,fac*set561_2_orig_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['561 to 647, no warping ' num2str(fac) 'x mag. residuals' ],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
fac = 2;
subplot(1,2,2)
hold on
plot(comb_set488_2_pos.x,comb_set488_2_pos.y,'g.')
quiver(comb_set488_2_pos.x,comb_set488_2_pos.y,fac*set488_2_orig_error_x,fac*set488_2_orig_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['488 to 647, no warping ' num2str(fac) 'x mag. residuals' ],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
saveas(gcf,[save_pth 'Visbead_errors_pre_warping.fig'],'fig')
saveas(gcf,[save_pth 'Visbead_errors_pre_warping.tif'],'tif')

% plot deviations after warping for right to left, stage to left registration
figure(2); clf
set(gcf,'position',[493 441 1212 482])
fac = 20;
subplot(1,2,1)
hold on
plot(comb_set561_2_pos.x,comb_set561_2_pos.y,'g.')
quiver(comb_set561_2_pos.x,comb_set561_2_pos.y,fac*set561_2_warp_error_x,fac*set561_2_warp_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['561 to 647, Order ' num2str(fac) ...
    'x mag. residuals'],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
fac = 20;
subplot(1,2,2)
hold on
plot(comb_set488_2_pos.x,comb_set488_2_pos.y,'g.')
quiver(comb_set488_2_pos.x,comb_set488_2_pos.y,fac*set488_2_warp_error_x,fac*set488_2_warp_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['488 to 647, Order' num2str(fac) ...
    'x mag. residuals'],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
saveas(gcf,[save_pth 'Visbead_errors_post_warping.fig'],'fig')
saveas(gcf,[save_pth 'Visbead_errors_post_warping.tif'],'tif')

% plot x and y offsets for right to left warping
figure(3); clf
set(gcf,'position',[493 441 1112 482])
subplot(1,2,1);
hist(set561_2_warp_error_x*nm_per_pix,15)
set(gca,'fontsize',ax_fs)
xlabel('measured x offsets [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title(['561 to 647 STD_x\times2.35=' num2str(std_set561_2_warp_error_x*2.35*nm_per_pix,'%0.1f') ' nm'],'fontsize',ti_fs)
set(gca,'position',[0.1544 0.2 0.3102 0.65])
subplot(1,2,2);
hist(set561_2_warp_error_y*nm_per_pix,15)
set(gca,'fontsize',ax_fs)
xlabel('measured y offsets [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title(['561 to 647 STD_y\times2.35=' num2str(std_set561_2_warp_error_y*2.35*nm_per_pix,'%0.1f') ' nm'],'fontsize',ti_fs)
set(gca,'position',[0.5948 0.2 0.3102 0.65])
saveas(gcf,[save_pth '561_2_647_xy_offsets.fig'],'fig')
saveas(gcf,[save_pth '561_2_647_xy_offsets.tif'],'tif')

% plot x and y offsets for right to left warping
figure(4); clf
set(gcf,'position',[493 441 1112 482])
subplot(1,2,1);
hist(set488_2_warp_error_x*nm_per_pix,15)
set(gca,'fontsize',ax_fs)
xlabel('measured x offsets [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title(['488 to 647 STD_x\times2.35=' num2str(std_set488_2_warp_error_x*2.35*nm_per_pix,'%0.1f') ' nm'],'fontsize',ti_fs)
set(gca,'position',[0.1544 0.2 0.3102 0.65])
subplot(1,2,2);
hist(set488_2_warp_error_y*nm_per_pix,15)
set(gca,'fontsize',ax_fs)
xlabel('measured y offsets [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title(['488 to 647 STD_y\times2.35=' num2str(std_set488_2_warp_error_y*2.35*nm_per_pix,'%0.1f') ' nm'],'fontsize',ti_fs)
set(gca,'position',[0.5948 0.2 0.3102 0.65])
saveas(gcf,[save_pth '488_2_647_xy_offsets.fig'],'fig')
saveas(gcf,[save_pth '488_2_647_xy_offsets.tif'],'tif')

% plot total offsets for left to right warping
figure(5); clf
set(gcf,'position',[493 441 1112 482])
subplot(1,2,1)
hist(nm_per_pix*sqrt(set561_2_warp_error_x.^2 + set561_2_warp_error_y.^2),40)
set(gca,'fontsize',ax_fs)
xlabel('measured displacements [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title('561to 647 measured displacements','fontsize',ti_fs)
subplot(1,2,2)
cdfplot(nm_per_pix*sqrt(set561_2_warp_error_x.^2 + set561_2_warp_error_y.^2))
set(gca,'fontsize',ax_fs)
xlabel('distance error [nm]','fontsize',la_fs)
ylabel('cumulative probability','fontsize',la_fs)
title('Cumulative distribution','fontsize',la_fs)
saveas(gcf,[save_pth '561_2_647_total_offsets.fig'],'fig')
saveas(gcf,[save_pth '561_2_647_total_offsets.tif'],'tif')

% plot total offsets for left to right warping
figure(6); clf
set(gcf,'position',[493 441 1112 482])
subplot(1,2,1)
hist(nm_per_pix*sqrt(set488_2_warp_error_x.^2 + set488_2_warp_error_y.^2),40)
set(gca,'fontsize',ax_fs)
xlabel('measured displacements [nm]','fontsize',la_fs)
ylabel('frequency','fontsize',la_fs)
title('488 to 647 measured displacements','fontsize',ti_fs)
subplot(1,2,2)
cdfplot(nm_per_pix*sqrt(set488_2_warp_error_x.^2 + set488_2_warp_error_y.^2))
set(gca,'fontsize',ax_fs)
xlabel('distance error [nm]','fontsize',la_fs)
ylabel('cumulative probability','fontsize',la_fs)
title('Cumulative distribution','fontsize',la_fs)
saveas(gcf,[save_pth '488_2_647_total_offsets.fig'],'fig')
saveas(gcf,[save_pth '488_2_647_total_offsets.tif'],'tif')
end
end