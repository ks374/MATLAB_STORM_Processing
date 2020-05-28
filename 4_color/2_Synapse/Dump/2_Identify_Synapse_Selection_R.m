%%
outpath = 'Z:\Chenghang\chenghaz_013_XB2_P4_Control_B\analysis\Result\';
load([outpath 'sizeshapematR_and_cent_water10_area_cutoff.mat']);
dsvoxel  =155/0.3;
%
centG = centGw;
sizeshape_matG = sizeshape_mat;
voxel = [15.5, 15.5, 70];
for i=1:numel(centG(:,1))
   centG(i,:) = centG(i,:).*voxel ;
end
rcentG = centG;

i=4;  % i = range(8)
sizeshapuse = sizeshape_mat;
val1w = (sizeshapuse(:,10+i)+1)./(sizeshapuse(:,2+i)+1);
val2w = log10((sizeshapuse(:,10+i)+1)*0.0158*0.0158*0.07);
numel(find(val1w<0.99))

Xn=80; Yn=60;
Xrange=[min(val1w) max(val1w)]; Yrange=[min(val2w) max(val2w)];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ; 
%
figure;
H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
close

cutoffg =60;
H2 = H;
H2(H2>cutoffg)=cutoffg;

figure;
pcolor(X,Y,H2)
%
savefig([outpath 'Rsizeshape_heatmap.fig'])
%%
load([outpath 'statsRwater10_area_cutoff.mat']);
%
figure;
pcolor(X,Y,H2)
% manually draw polygon on figure
currpoly=impoly
% return polygon coordinates
synapse_regiong=currpoly.getPosition
% save figure of selected polygon region
savefig([outpath 'synapse_selection_poly.fig'])

%Return centroid and stats lists for all clusters in selected polygon area
centGa2s=centG(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

centGa2ns=centG(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

rcentGa2s=rcentG(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

rcentGa2ns=rcentG(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

% return size/shape array for selected clusters
sizeshape_matGa2s=sizeshape_matG(find(inpolygon(val1w,val2w,...
    synapse_regiong(:,1),synapse_regiong(:,2))),:);

sizeshape_matGa2ns=sizeshape_matG(find(~inpolygon(val1w,val2w,...
    synapse_regiong(:,1),synapse_regiong(:,2))),:);
 
% return detailed pixel stats lists for selected clusters
statsGa2s=statsGwater(find(inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

statsGa2ns=statsGwater(find(~inpolygon(val1w,val2w,synapse_regiong(:,1),...
    synapse_regiong(:,2))),:);

size(centGa2s,1)

%
Real_non_select_area_thre = 4;
centGa2ns = centGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
rcentGa2ns = rcentGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
statsGa2ns = statsGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
sizeshape_matGa2ns = sizeshape_matGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);

%save([base_folder 'statslistG2sw10.mat'],'centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statslistG2nsw10.mat'],'centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
save([outpath 'statslistR2sw10.mat'],'centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statslistR2nsw10.mat'],'centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')

%
save([outpath 'statsR2sw10.mat'],'statsGa2s','centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statsR2nsw10.mat'],'statsGa2ns','centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')
%save([base_folder 'statsG2sw10.mat'],'statsGb2s','centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statsG2nsw10.mat'],'statsGb2ns','centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
