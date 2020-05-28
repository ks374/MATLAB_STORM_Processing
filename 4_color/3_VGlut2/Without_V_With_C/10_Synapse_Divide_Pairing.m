clear;clc
%%
base_folder = 'Z:\Chenghang\chenghaz_012_P8_Control_A\analysis\Result\5_V_Syn\';
outpath = base_folder;
C_path = 'Z:\Chenghang\chenghaz_012_P8_Control_A\analysis\Result\4_CTB\';
%%
load([C_path 'statslistC2sw10.mat']);
load([C_path 'statslistC2nsw10.mat']);
%
centRa2ns = centGa2ns;
centRa2s = centGa2s;
rcentRa2ns = rcentGa2ns;
rcentRa2s = rcentGa2s;
sizeshape_matRa2ns = sizeshape_matGa2ns;
sizeshape_matRa2s = sizeshape_matGa2s;
clear centGa2ns centGa2s rcentGa2ns rcentGa2s sizeshape_matGa2ns sizeshape_matGa2s
%%
load([base_folder 'statsSSSN2w10_edges_plus_C.mat']);
voxel = [15.5 15.5 70];

for i = 1:numel(statsRwater)
    rcentGa2s(i,:) = statsRwater(i).WeightedCentroid.*voxel;
    volumeGs(i) = statsRwater(i).Area;
end

centGs = rcentGa2s;%(rcentGa2s(:,2)<46000,:);
%

sizeshape_matRs = sizeshape_matRa2s;%(rcentPa2s(:,2)<46000,:);
centRs = rcentRa2s;%(rcentPa2s(:,2)<46000,:);

volumeRs = sizeshape_matRs(:,19);
%%
nn_Gs_Rs = zeros(size(centGs,1),1);

disp('startG2')
parfor i=1:size(centGs,1)
   nn_Gs_Rs(i) = min(pdist2(centGs(i,:),centRs));
%   nn_Gs_Pall(i) = min(pdist2(centGs(i,:),centP_all));
end
%
save([base_folder 'nearest_neightbor_pairing_gw10pw10.mat'],'nn_*')

%%
load([base_folder 'C_add_to_statsSSSNw10_edges_Vglut2.mat'],'tintsG_p140');

mints_g70s = (([tintsG_p140])./[volumeGs]);
%%
Xn=70; Yn=80; Xrange=[-0.2 3]; Yrange=[1.2 4.0];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ;
%
val1w = log10(mints_g70s +1)';
val2w = log10(nn_Gs_Rs);  
figure; H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange); close;
cutoffg = 20; 
H1 = H; H1(H1>cutoffg)=cutoffg;
figure; 
pcolor(X,Y,H1)

%%
k=200;
dataall = (cat(2,val1w,val2w));
datause = dataall(randi(numel(dataall(:,1)),[5000 1]),:);
[RDg,CDg,orderg]=optics(zscore(datause),k);
figure; plot(RDg(orderg))
% figure; plot(CDg(orderg))
%% get dbscan clusters from optics threshold and plot (geph)
Eps = 0.099; clustID = 1; classg = zeros(numel(RDg),1);
for i = 1:numel(RDg)
    if RDg(orderg(i))>Eps
        if CDg(orderg(i))<=Eps
            clustID = clustID + 1; classg(orderg(i)) = clustID;
        end
    else
        classg(orderg(i)) = clustID;
    end
end
max(classg)
figure; 
plot(datause(classg==3,1),datause(classg==3,2),'g.'); hold all
plot(datause(classg==2,1),datause(classg==2,2),'r.')
%%
figure; pcolor(X,Y,H1); hold on
%ezcontour(@(x,y)pdf(gm,[x y]),[0 180],[3 8]);
%sp(:,1) = (V(1,1)*val1wh + V(2,1)*val2wh)<=-40;
cls = ClassificationDiscriminant.fit(datause,classg);
%Plot the classification boundaries.
K = cls.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
L = cls.Coeffs(2,3).Linear;% boundary between the second and third classes                 
clear f
% Plot the curve K + [x,y]*L  = 0.
xval = [0 200];
yval = -((L(1)/L(2))*xval + K/L(2));
h2 = line(xval, yval);
set(h2,'Color','r','LineWidth',2)
savefig([base_folder 'storm_gs_ps_mints_nn_shell_2d_hist.fig'])

%%
pairedg_idx = -((L(1)/L(2))*val1w + K/L(2))>val2w;
numel(find(pairedg_idx))/numel((pairedg_idx))
% pairedp_idx = -((L2(1)/L2(2))*val5w + K2/L2(2))>val6w;
% numel(find(pairedp_idx))/numel((pairedp_idx))
%
save([base_folder 'pairing_index_ps_gs_withedges.mat'],'pairedg_idx')
% save([base_folder 'wga_normalization_pixels_withedges.mat'],'CW2')

%%
load([outpath 'statsSSSN2w10_edges_plus_C.mat']);
statsRwater_sssns = statsRwater(pairedg_idx);
statsRwater_sssnn = statsRwater(~pairedg_idx);
save([base_folder 'R_paired_C.mat'],'statsRwater_sssns','statsRwater_sssnn');
