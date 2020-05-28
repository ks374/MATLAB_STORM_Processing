%function f = gen_bead_warp(arg1)
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
%arg1 = ['Y:/backup2/ysigal/Brad_samples/' ...
%    '102514_P40G_750_MBP_647_Phal_488_WGA/coverslips/' ...
%    '102514_P40G_750_MBP_647_Phal_488_WGA_cs01/'];
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

% file i/o stuff
if exist(save_pth,'dir') == 0
    mkdir(save_pth)
end
numIRbeadmov = numel(dir([acq_path 'IRbeads*' '*.dax']));
numVisbeadmov = numel(dir([acq_path 'Visbeads*' '*.dax']));
if numIRbeadmov>1
IRworkspace_nm_r2 = [save_pth 'IRbeadworkspace_region2_' num2str(scale) 'x.mat'];
track_file_nms = dir([base_pth 'IRbeads*' '_2_' '*mlist.bin'] );
if isempty(track_file_nms)
    track_file_nms = dir([base_pth 'IRbeads*' '*mlist.bin'] );
end
if exist(IRworkspace_nm_r2,'file') == 0 
for p = 1:5
    if p == 1
        len = 1;
        match_radius = 12*scale; % radius tolerance for matching left/right channels

    elseif p==2
        len = 2;
        match_radius = 7*scale; % radius tolerance for matching left/right channels
    elseif p==3
        len = 5;
        match_radius = 3*scale; % radius tolerance for matching left/right channels
    elseif p==4
        len = 10;
        match_radius = 1*scale; % radius tolerance for matching left/right channels

    else
        len = length(track_file_nms);
        match_radius = 0.4*scale; % radius tolerance for matching left/right channels

    end
        

for k = 1:len
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

    if p==1
         tform_750start = maketform('affine',[1 0 0; 0 1 0; -init_750xshift -init_750yshift 1]);
    else
        tform_750start = tform_750_2_647;
    end

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
if p==1
     tform_750_2_647 = cp2tform((set750_1_points),(set750_2_points),'similarity');
else

try
    polyorder_R2L = 4;
    tform_750_2_647 = cp2tform(set750_1_points,set750_2_points,'polynomial',polyorder_R2L); % why sometimes error??
    catch
        try
        polyorder_R2L = 3;
        tform_750_2_647 = cp2tform(set750_1_points,set750_2_points,'polynomial',polyorder_R2L); % why sometimes error??
        catch
        polyorder_R2L = 2;
        tform_750_2_647 = cp2tform(set750_1_points,set750_2_points,'polynomial',polyorder_R2L); % why sometimes error??
        try
        tform_750_2_647 = cp2tform((set750_1_points),(set750_2_points),'similarity'); % why sometimes error??
    
        end
    end
end
end
[warped_set750_2_pos.x,warped_set750_2_pos.y] = tforminv(tform_750_2_647,comb_set750_2_pos.x,comb_set750_2_pos.y);

set750_2_warp_error_x = comb_set750_1_pos.x - warped_set750_2_pos.x;
set750_2_warp_error_y = comb_set750_1_pos.y - warped_set750_2_pos.y;
std_set750_2_warp_error_x = std(set750_2_warp_error_x);
std_set750_2_warp_error_y = std(set750_2_warp_error_y);
std_total750_warp_error = std(sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2));


set750_2_orig_error_x = comb_set750_1_pos.x - comb_set750_2_pos.x - init_750xshift;
set750_2_orig_error_y = comb_set750_1_pos.y - comb_set750_2_pos.y - init_750yshift;
std_set750_2_orig_error_x = std(set750_2_orig_error_x);
std_set750_2_orig_error_y = std(set750_2_orig_error_y);

[cdf750(p).y, cdf750(p).x] = ecdf(nm_per_pix*sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2));
cdf90_750 = (cdf750(p).x(find(cdf750(p).y>0.9,1,'first')));
disp(['90% of 750 beads aligned to ' num2str(cdf90_750) 'nm ,using ' num2str(length(set750_1_points(:,1))) ' beads'])
if p>=2
   if (cdf750(p).x(find(cdf750(p).y>0.9,1,'first')))== (cdf750(p-1).x(find(cdf750(p-1).y>0.9,1,'first')))
   disp('done!')
   break
   end
end

end

save(IRworkspace_nm_r2,'tform_750_2_647');
end
end
% file i/o stuff
if numVisbeadmov>1

Vis_workspace_nm_r2 = [save_pth 'Visbeadworkspace_region2_' num2str(scale) 'x.mat'];
track_file_nms = dir([base_pth 'Visbeads*' '_2_' '*mlist.bin'] );
if isempty(track_file_nms)
    track_file_nms = dir([base_pth 'Visbeads*' '*mlist.bin'] );
end
if exist(Vis_workspace_nm_r2,'file') == 0 
for p = 1:5
    if p == 1
        len = 1;
        match_radius = 12*scale; % radius tolerance for matching left/right channels

    elseif p==2
        len = 2;
        match_radius = 7*scale; % radius tolerance for matching left/right channels
    elseif p==3
        len = 5;
        match_radius = 2*scale; % radius tolerance for matching left/right channels
    elseif p==4
        len = 10;
        match_radius = 1*scale; % radius tolerance for matching left/right channels

    else
        len = length(track_file_nms);
        match_radius = 0.5*scale; % radius tolerance for matching left/right channels

    end

for k = 1:len
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

    if p==1
         tform_561start = maketform('affine',[1 0 0; 0 1 0; -init_561xshift -init_561yshift 1]);
    else
        tform_561start = tform_561_2_647;
    end
    
    if p==1
         tform_488start = maketform('affine',[1 0 0; 0 1 0; -init_488xshift -init_488yshift 1]);
    else
        tform_488start = tform_488_2_647;
    end

 [A(k).matched561,A(k).unmatched561] = corr_mols(A(k).set561_pos,A(k).set647_pos,tform_561start,match_radius);
 [aa ab ac] = unique(A(k).matched561.set2_inds);
 A(k).matched561.set1_inds = A(k).matched561.set1_inds(ab);
 A(k).matched561.set2_inds = A(k).matched561.set2_inds(ab);

 %   disp([num2str(2*length(A(k).matched561.set1_inds)) '/' ...
 %       num2str(2*length(A(k).matched561.set1_inds) ...
 %      +length(A(k).unmatched561.set1_inds)+ ...
 %       length(A(k).unmatched561.set2_inds)) ' 561 molecules matched'])
    
  [A(k).matched488,A(k).unmatched488] = corr_mols(A(k).set488_pos,A(k).set647_pos,tform_488start,match_radius);
 [aa ab ac] = unique(A(k).matched488.set2_inds);
 A(k).matched488.set1_inds = A(k).matched488.set1_inds(ab);
 A(k).matched488.set2_inds = A(k).matched488.set2_inds(ab);
    
 %   disp([num2str(2*length(A(k).matched488.set1_inds)) '/' ...
 %       num2str(2*length(A(k).matched488.set1_inds) ...
 %       +length(A(k).unmatched488.set1_inds)+ ...
 %       length(A(k).unmatched488.set2_inds)) ' 488 molecules matched'])
 %  disp(['bead image number' num2str(k)])
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

if p==1
     tform_561_2_647 = cp2tform((set561_1_points),(set561_2_points),'similarity');
else

try
    polyorder_R2L = 4;
    tform_561_2_647 = cp2tform(set561_1_points,set561_2_points,'polynomial',polyorder_R2L); % why sometimes error??
catch
    try
        polyorder_R2L = 3;
        tform_561_2_647 = cp2tform(set561_1_points,set561_2_points,'polynomial',polyorder_R2L); % why sometimes error??
    catch
        polyorder_R2L = 2;
        tform_561_2_647 = cp2tform(set561_1_points,set561_2_points,'polynomial',polyorder_R2L); % why sometimes error??
    end
end
end

if p==1
     tform_488_2_647 = cp2tform((set488_1_points),(set488_2_points),'similarity');
else

try
    polyorder_R2L = 4;
    tform_488_2_647 = cp2tform(set488_1_points,set488_2_points,'polynomial',polyorder_R2L); % why sometimes error??
catch
    try
        polyorder_R2L = 3;
        tform_488_2_647 = cp2tform(set488_1_points,set488_2_points,'polynomial',polyorder_R2L); % why sometimes error??
    catch
        polyorder_R2L = 2;
        tform_488_2_647 = cp2tform(set488_1_points,set488_2_points,'polynomial',polyorder_R2L); % why sometimes error??
    end
end
end


[warped_set561_2_pos.x,warped_set561_2_pos.y] = tforminv(tform_561_2_647,comb_set561_2_pos.x,comb_set561_2_pos.y);

set561_2_warp_error_x = comb_set561_1_pos.x - warped_set561_2_pos.x;
set561_2_warp_error_y = comb_set561_1_pos.y - warped_set561_2_pos.y;
std_set561_2_warp_error_x = std(set561_2_warp_error_x);
std_set561_2_warp_error_y = std(set561_2_warp_error_y);
std_total561_warp_error = std(sqrt(set561_2_warp_error_x.^2 + set561_2_warp_error_y.^2));


set561_2_orig_error_x = comb_set561_1_pos.x - comb_set561_2_pos.x - init_561xshift;
set561_2_orig_error_y = comb_set561_1_pos.y - comb_set561_2_pos.y - init_561yshift;
std_set561_2_orig_error_x = std(set561_2_orig_error_x);
std_set561_2_orig_error_y = std(set561_2_orig_error_y);


[warped_set488_2_pos.x,warped_set488_2_pos.y] = tforminv(tform_488_2_647,comb_set488_2_pos.x,comb_set488_2_pos.y);

set488_2_warp_error_x = comb_set488_1_pos.x - warped_set488_2_pos.x;
set488_2_warp_error_y = comb_set488_1_pos.y - warped_set488_2_pos.y;
std_set488_2_warp_error_x = std(set488_2_warp_error_x);
std_set488_2_warp_error_y = std(set488_2_warp_error_y);
std_total488_warp_error = std(sqrt(set488_2_warp_error_x.^2 + set488_2_warp_error_y.^2));

set488_2_orig_error_x = comb_set488_1_pos.x - comb_set488_2_pos.x - init_488xshift;
set488_2_orig_error_y = comb_set488_1_pos.y - comb_set488_2_pos.y - init_488yshift;
std_set488_2_orig_error_x = std(set488_2_orig_error_x);
std_set488_2_orig_error_y = std(set488_2_orig_error_y);

[cdf488(p).y, cdf488(p).x] = ecdf(nm_per_pix*sqrt(set488_2_warp_error_x.^2 + set488_2_warp_error_y.^2));
cdf90_488 = (cdf488(p).x(find(cdf488(p).y>0.9,1,'first')));
disp(['90% of 488 beads aligned to ' num2str(cdf90_488) 'nm ,using ' num2str(length(set488_1_points(:,1))) ' beads'])

[cdf561(p).y, cdf561(p).x] = ecdf(nm_per_pix*sqrt(set561_2_warp_error_x.^2 + set561_2_warp_error_y.^2));
cdf90_561 = (cdf561(p).x(find(cdf561(p).y>0.9,1,'first')));
disp(['90% of 561 beads aligned to ' num2str(cdf90_561) 'nm ,using ' num2str(length(set488_1_points(:,1))) ' beads'])
if p>=2
   if (cdf561(p).x(find(cdf561(p).y>0.9,1,'first')))== (cdf561(p-1).x(find(cdf561(p-1).y>0.9,1,'first'))) & ...
        (cdf488(p).x(find(cdf488(p).y>0.9,1,'first')))==(cdf488(p-1).x(find(cdf488(p-1).y>0.9,1,'first')))
   disp('done!')
   break
   end
   end
end
%
save(Vis_workspace_nm_r2,'tform_561_2_647','tform_488_2_647');
end
end
%save output plots to demonstrate goodness of fit
if exist([save_pth 'Cumulative_distribution_for_registration_self.tif'],'file') == 0
%%plot deviations before warping for right to left, stage to left registration
figure(1); clf
set(gcf,'position',[493 441 1512 882])
fac = 2;
subplot(2,3,2)
hold on
if numVisbeadmov>1
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
subplot(2,3,1)
hold on
plot(comb_set488_2_pos.x,comb_set488_2_pos.y,'g.')
quiver(comb_set488_2_pos.x,comb_set488_2_pos.y,fac*set488_2_orig_error_x,fac*set488_2_orig_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['488 to 647, no warping ' num2str(fac) 'x mag. residuals' ],'fontsize',ti_fs)
%self_warping
fac = 20;
subplot(2,3,5)
hold on
plot(comb_set561_2_pos.x,comb_set561_2_pos.y,'g.')
quiver(comb_set561_2_pos.x,comb_set561_2_pos.y,fac*set561_2_warp_error_x,fac*set561_2_warp_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['561 to 647, Order ' num2str(polyorder_R2L) ' poly. warp (self), ' num2str(fac) ...
    'x mag. residuals'],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
fac = 20;
subplot(2,3,4)
hold on
plot(comb_set488_2_pos.x,comb_set488_2_pos.y,'g.')
quiver(comb_set488_2_pos.x,comb_set488_2_pos.y,fac*set488_2_warp_error_x,fac*set488_2_warp_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['488 to 647, Order ' num2str(polyorder_R2L) ' poly. warp (self), ' num2str(fac) ...
    'x mag. residuals'],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
end
if numIRbeadmov>1
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
fac = 2;
subplot(2,3,3)
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

fac = 20;
subplot(2,3,6)
hold on
plot(comb_set750_2_pos.x,comb_set750_2_pos.y,'g.')
quiver(comb_set750_2_pos.x,comb_set750_2_pos.y,fac*set750_2_warp_error_x,fac*set750_2_warp_error_y,0,'r');
axis ij
set(gca,'fontsize',ax_fs)
xlabel('x axis [pixel]','fontsize',la_fs)
ylabel('y axis [pixel]','fontsize',la_fs)
title(['750 to 647, Order ' num2str(polyorder_R2L) ' poly. warp (self), ' num2str(fac) ...
    'x mag. residuals'],'fontsize',ti_fs)
hold off
xlim([-10 265*scale]); ylim([-20 265*scale])
end
saveas(gcf,[save_pth '_BeadMap_for_registration_self.tif'],'tif')
saveas(gcf,[save_pth '_BeadMap_for_registration_self.fig'],'fig')

% plot total offsets for left to right warping
figure(2); clf
set(gcf,'position',[493 441 1312 1082])
if numVisbeadmov>1

subplot(1,3,2)
cdfplot(nm_per_pix*sqrt(set561_2_warp_error_x.^2 + set561_2_warp_error_y.^2))
set(gca,'fontsize',ax_fs)
xlabel('distance error [nm]','fontsize',la_fs)
ylabel('cumulative probability','fontsize',la_fs)
title({'Cumulative distribution self-warping 561 to 647'; ['STD=' num2str(std_total561_warp_error*nm_per_pix,'%0.1f') 'nm']; ['90% of 561 beads aligned to ' num2str(cdf90_561) 'nm']} ,'fontsize',la_fs)
subplot(1,3,1)
cdfplot(nm_per_pix*sqrt(set488_2_warp_error_x.^2 + set488_2_warp_error_y.^2))
set(gca,'fontsize',ax_fs)
xlabel('distance error [nm]','fontsize',la_fs)
ylabel('cumulative probability','fontsize',la_fs)
title({'Cumulative distribution self-warping 488 to 647'; ['STD=' num2str(std_total488_warp_error*nm_per_pix,'%0.1f') 'nm']; ['90% of 488 beads aligned to ' num2str(cdf90_488) 'nm']} ,'fontsize',la_fs)
end
if numIRbeadmov>1

subplot(1,3,3)
cdfplot(nm_per_pix*sqrt(set750_2_warp_error_x.^2 + set750_2_warp_error_y.^2))
set(gca,'fontsize',ax_fs)
xlabel('distance error [nm]','fontsize',la_fs)
ylabel('cumulative probability','fontsize',la_fs)
title({'Cumulative distribution self-warping 750 to 647'; ['STD=' num2str(std_total750_warp_error*nm_per_pix,'%0.1f') 'nm']; ['90% of 750 beads aligned to ' num2str(cdf90_750) 'nm']} ,'fontsize',la_fs)
end
saveas(gcf,[save_pth 'Cumulative_distribution_for_registration_self.fig'],'fig')
saveas(gcf,[save_pth 'Cumulative_distribution_for_registration_self.tif'],'tif')
end


