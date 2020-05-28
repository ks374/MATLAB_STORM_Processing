clear all

%first separate synaptic groups,for this pass, 
%first  do it quickly with manual separation of synaptic and non_synaptic
%%
%base_folder = 'Y:/backup/ysigal/GLD3_bistratified_RGC/elastic_thresh/';
base_folder = '/n/contefs1/backup/ysigal/GLD3_bistratified_RGC/elastic_thresh/';
base_folder1 = '/n/contefs1/backup/ysigal/GLD3_bistratified_RGC/elastic_align/';
%load in presynaptic cluster list
%load([base_folder 'sizeshapematP_and_cent_log_t2w30.mat']);
%
load([base_folder 'statslistP2sw10.mat']);
load([base_folder 'statslistP2nsw10.mat']);
load([base_folder 'statslistG2sw10.mat']);
load([base_folder 'statslistG2nsw10.mat']);

sizeshape_matGns = sizeshape_matGa2ns;%(rcentGa2ns(:,2)<46000,:);
centGa2ns2 = centGa2ns;%(rcentGa2ns(:,2)<46000,:);
centGns = rcentGa2ns;%(rcentGa2ns(:,2)<46000,:);

sizeshape_matGs = sizeshape_matGa2s;%(rcentGa2s(:,2)<46000,:);
centGa2s2 = centGa2s;%(rcentGa2s(:,2)<46000,:);
centGs = rcentGa2s;%(rcentGa2s(:,2)<46000,:);
%
sizeshape_matPns = sizeshape_matPa2ns;%(rcentPa2ns(:,2)<46000,:);
centPa2ns2 = centPa2ns;%(rcentPa2ns(:,2)<46000,:);
centPns = rcentPa2ns;%(rcentPa2ns(:,2)<46000,:);

sizeshape_matPs = sizeshape_matPa2s;%(rcentPa2s(:,2)<46000,:);
centPa2s2 = centPa2s;%(rcentPa2s(:,2)<46000,:);
centPs = rcentPa2s;%(rcentPa2s(:,2)<46000,:);

centP_all = cat(1, centPs, centPns);
volumePs = sizeshape_matPs(:,19);
centG_all = cat(1, centGns, centGs);
volumeGs = sizeshape_matGs(:,19);

%
load([base_folder1 'rotated_ds_data.mat']);
dsvoxel  =158/0.3;
clear CW
BW3 = BP3;
[CW(:,1), CW(:,2), CW(:,3)] = ind2sub(size(BW3),find(BW3));
for j=1:3
    CW(:,j) = CW(:,j).*dsvoxel;
end
regionbounds = regionbounds_ipl;
CW = CW(CW(:,1)>(regionbounds(1,1)*dsvoxel(1)) & ...
    CW(:,1)<(regionbounds(2,1)*dsvoxel(1)) & ...
    CW(:,2)>(regionbounds(1,2)*dsvoxel(1)) & ...
    CW(:,2)<(regionbounds(2,2)*dsvoxel(1)) & ...
    CW(:,3)>(regionbounds(1,3)*dsvoxel(1)) & ...
    CW(:,3)<(regionbounds(2,3)*dsvoxel(1)),:);
CW(:,1)=CW(:,1)-regionbounds(1,1)*dsvoxel;
CW(:,2)=CW(:,2)-regionbounds(1,2)*dsvoxel;
CW(:,3)=CW(:,3)-regionbounds(1,3)*dsvoxel;
%
kk = 86;
lamsize = max(CW(:,1))-min(CW(:,1));
samp_freq = kk;
binedges = 0:lamsize/samp_freq:lamsize;
hw1 = histc(CW(:,1),binedges);%
CW2 = CW;%(CW(:,1)>11000,:);
%CW2 = CW2(CW2(:,1)<46000,:);
%
centGs_rand = CW2(randi(numel(CW2(:,1)),(numel(centGs(:,1))),1),:);
centGs_rand = centGs_rand(:,[2 1 3]);
centGs_rand = centGs_rand + randi(round(dsvoxel),numel(centGs(:,1)),3); 
centGall_rand = CW2(randi(numel(CW2(:,1)),numel(centG_all(:,1)),1),:);
centGall_rand = centGall_rand(:,[2 1 3]);
centGall_rand = centGall_rand + randi(round(dsvoxel),numel(centG_all(:,1)),3); 

centPs_rand = CW2(randi(numel(CW2(:,1)),round(numel(centPs(:,1))),1),:);
centPs_rand = centPs_rand(:,[2 1 3]);
centPs_rand = centPs_rand + randi(round(dsvoxel),numel(centPs(:,1)),3); 

centPall_rand = CW2(randi(numel(CW2(:,1)),numel(centP_all(:,1)),1),:);
centPall_rand = centPall_rand(:,[2 1 3]);
centPall_rand = centPall_rand + randi(round(dsvoxel),numel(centP_all(:,1)),3); 

%
nn_Gall_Ps = zeros(size(centG_all,1),1);
nn_Garand_Psrand = zeros(size(centG_all,1),1);
nn_Gs_Ps = zeros(size(centGs,1),1);
nn_Gns_Ps = zeros(size(centGns,1),1);
%nn_Gall_Pall = zeros(size(centG_all,1),1);
%nn_Grand_Pall = zeros(size(centG_all,1),1);
%nn_Gs_Pall = zeros(size(centGs,1),1);
%nn_Gns_Pall = zeros(size(centGns,1),1);

nn_Pall_Gs = zeros(size(centP_all,1),1);
nn_Prand_Gsrnad = zeros(size(centP_all,1),1);
nn_Ps_Gs = zeros(size(centPs,1),1);
nn_Pns_Gs = zeros(size(centPns,1),1);
%nn_Pall_Gall = zeros(size(centP_all,1),1);
%nn_Prand_Gall = zeros(size(centP_all,1),1);
%nn_Ps_Gall = zeros(size(centPs,1),1);
%nn_Pns_Gall = zeros(size(centPns,1),1);
%
parpool(32)
disp('startG1')
parfor i=1:size(centG_all,1)
nn_Gall_Ps(i) = min(pdist2(centG_all(i,:),centPs));
nn_Grand_Ps(i) = min(pdist2(centGall_rand(i,:),centPs_rand));
%nn_Gall_Pall(i) = min(pdist2(centG_all(i,:),centP_all));
%nn_Grand_Pall(i) = min(pdist2(centGall_rand(i,:),centP_all));
end
%
disp('startG2')
parfor i=1:size(centGs,1)
   nn_Gs_Ps(i) = min(pdist2(centGs(i,:),centPs));
%   nn_Gs_Pall(i) = min(pdist2(centGs(i,:),centP_all));
end
disp('startG3')
parfor i=1:size(centGns,1)
   nn_Gns_Ps(i) = min(pdist2(centGns(i,:),centPs));
%   nn_Gns_Pall(i) = min(pdist2(centGns(i,:),centP_all));
end
%
disp('startP1')
parfor i=1:size(centP_all,1)
nn_Pall_Gs(i) = min(pdist2(centP_all(i,:),centGs));
nn_Prand_Gs(i) = min(pdist2(centPall_rand(i,:),centGs_rand));
%nn_Pall_Gall(i) = min(pdist2(centP_all(i,:),centG_all));
%nn_Prand_Gall(i) = min(pdist2(centPall_rand(i,:),centG_all));
end
%
disp('startP2')
parfor i=1:size(centPs,1)
   nn_Ps_Gs(i) = min(pdist2(centPs(i,:),centGs));
 %  nn_Ps_Gall(i) = min(pdist2(centPs(i,:),centG_all));
end
disp('startP3')
parfor i=1:size(centPns,1)
   nn_Pns_Gs(i) = min(pdist2(centPns(i,:),centGs));
 %  nn_Pns_Gall(i) = min(pdist2(centPns(i,:),centG_all));
end
%
save([base_folder 'nearest_neightbor_pairing_gw10pw10.mat'],'nn_*')
%
binsl = 0:0.05:5;
[hy1, hx2] = hist(log10(nn_Gall_Ps),binsl);
[hy2, hx2] = hist(log10(nn_Gs_Ps),binsl);
[hy3, hx2] = hist(log10(nn_Gns_Ps),binsl);
[hy4, hx2] = hist(log10(nn_Grand_Ps),binsl);
%[hy1a, hx2] = hist(log10(nn_Gall_Pall),binsl);
%[hy2a, hx2] = hist(log10(nn_Gs_Pall),binsl);
%[hy3a, hx2] = hist(log10(nn_Gns_Pall),binsl);
%[hy4a, hx2] = hist(log10(nn_Grand_Pall),binsl);

[hy5, hx2] = hist(log10(nn_Pall_Gs),binsl);
[hy6, hx2] = hist(log10(nn_Ps_Gs),binsl);
[hy7, hx2] = hist(log10(nn_Pns_Gs),binsl);
[hy8, hx2] = hist(log10(nn_Prand_Gs),binsl);
%[hy5a, hx2] = hist(log10(nn_Pall_Gall),binsl);
%[hy6a, hx2] = hist(log10(nn_Ps_Gall),binsl);
%[hy7a, hx2] = hist(log10(nn_Pns_Gall),binsl);
%[hy8a, hx2] = hist(log10(nn_Prand_Gall),binsl);
%
figure;
plot(hx2,hy1,'k'); alpha(0.5); hold on
plot(hx2,hy2,'g'); alpha(0.5); hold on
plot(hx2,hy3,'r'); alpha(0.5); hold on
plot(hx2,hy4*0.5,'b'); alpha(0.5)
savefig([base_folder 'nnGsplit_Ps_counts_with_rand_and_all.fig'])
%
%figure;
%plot(hx2,hy1a,'k'); alpha(0.5); hold on
%plot(hx2,hy2a,'g'); alpha(0.5); hold on
%plot(hx2,hy3a,'r'); alpha(0.5); hold on
%plot(hx2,hy4a*0.4,'b'); alpha(0.5)
%savefig([base_folder 'nnGsplit_Pall_counts_with_rand_and_all.fig'])

figure;
plot(hx2,hy5,'k'); alpha(0.5); hold on
plot(hx2,hy6,'g'); alpha(0.5); hold on
plot(hx2,hy7,'r'); alpha(0.5); hold on
plot(hx2,hy8*0.5,'b'); alpha(0.5)
savefig([base_folder 'nnPsplit_Gs_counts_with_rand_and_all.fig'])
%
%figure;
%plot(hx2,hy5a,'k'); alpha(0.5); hold on
%plot(hx2,hy6a,'g'); alpha(0.5); hold on
%plot(hx2,hy7a,'r'); alpha(0.5); hold on
%plot(hx2,hy8a*0.4,'b'); alpha(0.5)
%savefig([base_folder 'nnPsplit_Gall_counts_with_rand_and_all.fig'])

%%
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%
clear hyall_use hys_use hyns_use
hyall_use = cat(1,hy1,hy5);
hys_use = cat(1,hy2,hy6);
hyns_use = cat(1,hy3,hy7);
for i=1:2
    disp(i)
[xData, yData] = prepareCurveData( hx2, hyall_use(i,:));
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts. Lower = [0 1.8 0 0 2.4 0];
opts.StartPoint = [12849 2.2 0.3 i^2*15803.883942715 2.8 0.4];
opts. Upper = [Inf 2.3 0.4 Inf 3.5 0.5];

[fr, gof] = fit( xData, yData, ft, opts );
confintall = confint(fr);
%
%figure( 'Name', 'untitled fit 1' );
%subplot(1,3,1)
%h = plot( fr, xData, yData );
clear xData yData
[xData, yData] = prepareCurveData( hx2, hys_use(i,:));
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts. Lower = [0 fr.b1 fr.c1 0 fr.b2 fr.c2];
opts.StartPoint = [hys_use(i,round(fr.b1*10)) fr.b1 fr.c1 hys_use(i,round(fr.b2*10)) fr.b2 fr.c2];
opts. Upper = [Inf fr.b1*1.001 fr.c1*1.001 Inf fr.b2*1.001 fr.c2*1.001];
[fitresults, gofs] = fit( xData, yData, ft, opts );
%subplot(1,3,2)
%h = plot( fitresults, xData, yData );


[xData2, yData2] = prepareCurveData( hx2, hyns_use(i,:));
clear ft opts
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts. Lower = [0 fr.b1 fr.c1 0 fr.b2 fr.c2];
opts.StartPoint = [hyns_use(i,round(fr.b1*10)) fr.b1 fr.c1 hyns_use(i,round(fr.b2*10)) fr.b2 fr.c2];
opts. Upper = [Inf fr.b1*1.001 fr.c1*1.001 Inf fr.b2*1.001 fr.c2*1.001];
[fitresultns, gofns] = fit( xData2, yData2, ft, opts );
%subplot(1,3,3)
%h = plot( fitresultns, xData2, yData2 );
goodfit(i,1) = gof.adjrsquare;
goodfit(i,2) = gofs.adjrsquare;
goodfit(i,3) = gofns.adjrsquare;
%
sa1(i) = fitresults.a1;
sb1(i) = fitresults.b1;
sc1(i) = fitresults.c1;
sa2(i) = fitresults.a2;
sb2(i) = fitresults.b2;
sc2(i) = fitresults.c2;
nsa1(i) = fitresultns.a1;
nsb1(i) = fitresultns.b1;
nsc1(i) = fitresultns.c1;
nsa2(i) = fitresultns.a2;
nsb2(i) = fitresultns.b2;
nsc2(i) = fitresultns.c2;
syngeph(i) =  fitresults.a1*fitresults.c1 *10*sqrt(2*pi);
paired_s(i) = fitresults.a1*fitresults.c1 *10*sqrt(2*pi)/...
    ( fitresults.a1*fitresults.c1 *10*sqrt(2*pi)+ ...
    fitresults.a2*fitresults.c2 *10*sqrt(2*pi));
unpaired_s(i) = fitresults.a2*fitresults.c2 *10*sqrt(2*pi)/...
    ( fitresults.a1*fitresults.c1 *10*sqrt(2*pi)+ ...
    fitresults.a2*fitresults.c2 *10*sqrt(2*pi));
paired_ns(i) = fitresultns.a1*fitresultns.c1 *10*sqrt(2*pi)/...
    ( fitresultns.a1*fitresultns.c1 *10*sqrt(2*pi)+ ...
    fitresultns.a2*fitresultns.c2 *10*sqrt(2*pi));
unpaired_ns(i) = fitresultns.a2*fitresultns.c2 *10*sqrt(2*pi)/...
    ( fitresultns.a1*fitresultns.c1 *10*sqrt(2*pi)+ ...
    fitresultns.a2*fitresultns.c2 *10*sqrt(2*pi));
end
%%
figure;
NumTrials=5; % Number of trail fits per sample
PeakShape=1; % Gaussian=1, Lorentzian=2, Logistic=3, Pearson=4
% "extra" determines the shape of the Pearson function. When extra=1
% the shape is Lorentzian; when extra is large (>20), the shape 
% approaches Gaussian. At small values (<1), the shape is a pointy "cusp".
extra=2; 
center=5; % x-value in center of fitted range
window =2000;
NumPeaks = 2;
%startvector=[100 50 300 100 ];
startvector=[2 1 2.7 1];

[FitResultsG,LowestError]=peakfit([hx2' (hy5)'],center,window,NumPeaks,PeakShape,extra,NumTrials,startvector); 

%%
goodfit
paired_s
unpaired_s
paired_ns
unpaired_ns
%%
load([base_folder 'add_to_statsGw10_edges.mat'],'tintsG_p140','tintsG_p220');
tintsG_p1402 = tintsG_p140;%(rcentGa2s(:,2)<46000);
tintsG_p2202 = tintsG_p220;%(rcentGa2s(:,2)<46000);
x = load([base_folder 'add_to_statsP2w10.mat'],'tintsG_p70');
tintsP_g70 = x.tintsG_p70;
tintsP_g702 = tintsP_g70%(rcentPa2s(:,2)<46000);

mints_g140s = (([tintsG_p1402])./[volumeGs]');
mints_g220s = (([tintsG_p2202])./[volumeGs]');
mints_p70s = (([tintsP_g702])./[volumePs]');
%%
Xn=70; Yn=80; Xrange=[-0.2 3]; Yrange=[1.2 3.2];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ;
%
val1w = log10(mints_g220s +1)';
val2w = log10(nn_Gs_Ps);  
figure; H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange); close;
cutoffg =500; H1 = H; H1(H1>cutoffg)=cutoffg;
figure; pcolor(X,Y,H1)

%%
val5w = log10(mints_p70s +1)';
val6w = log10(nn_Ps_Gs);  
figure; H = hist2d(cat(2,val5w,val6w),Xn,Yn,[-0.2 2],Yrange); close;
cutoffg =500; H3 = H; H3(H3>cutoffg)=cutoffg;
figure; pcolor(X,Y,H3)
%savefig([path 'storm_ps_gs_mints_nn_shell_2d_hist.fig'])
%%
k=200;
dataall = (cat(2,val1w,val2w));
datause = dataall(randi(numel(dataall(:,1)),[5000 1]),:);
[RDg,CDg,orderg]=optics(zscore(datause),k);
figure; plot(RDg(orderg))
%% get dbscan clusters from optics threshold and plot (geph)
Eps = 0.69; clustID = 1; classg = zeros(numel(RDg),1);
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
k=200;
dataall = (cat(2,val5w,val6w));
datause = dataall(randi(numel(dataall(:,1)),[5000 1]),:);
[RDg,CDg,orderg]=optics(zscore(datause),k);
figure; plot(RDg(orderg))
%% get dbscan clusters from optics threshold and plot (geph)
Eps = 0.25; clustID = 1; classg = zeros(numel(RDg),1);
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
figure; pcolor(X,Y,H3); hold on
%ezcontour(@(x,y)pdf(gm,[x y]),[0 180],[3 8]);
%sp(:,1) = (V(1,1)*val1wh + V(2,1)*val2wh)<=-40;
cls = ClassificationDiscriminant.fit(datause,classg);
%Plot the classification boundaries.
K2 = cls.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
L2 = cls.Coeffs(2,3).Linear;% boundary between the second and third classes                 
clear f
% Plot the curve K + [x,y]*L  = 0.
xval = [0 200];
yval = -((L2(1)/L2(2))*xval + K2/L2(2));
h2 = line(xval, yval);
set(h2,'Color','r','LineWidth',2)
savefig([base_folder 'storm_ps_gs_mints_nn_shell_2d_hist.fig'])
%%
pairedg_idx = -((L(1)/L(2))*val1w + K/L(2))>val2w;
numel(find(pairedg_idx))/numel((pairedg_idx))
pairedp_idx = -((L2(1)/L2(2))*val5w + K2/L2(2))>val6w;
numel(find(pairedp_idx))/numel((pairedp_idx))
%%
save([base_folder 'pairing_index_ps_gs_withedges.mat'],'pairedg_idx','pairedp_idx')
save([base_folder 'wga_normalization_pixels_withedges.mat'],'CW2')
%%
max(binedges)
%%
figure;
currcent = centGs(pairedg_idx,:);
currstat = sizeshape_matGs(pairedg_idx,:);

volume = currstat(:,19);
tints = currstat(:,20);
clear zsize  ztints zdens numbins hy2

for i=1:(numel(binedges))
    if i==(numel(binedges))
        zsize(i) = mean(volume(currcent(:,2)== binedges(i)));
        zdens(i) = numel([currcent(currcent(:,2)== binedges(i),2)]);
        ztints(i) = sum(tints(currcent(:,2)== binedges(i)));
    else
        zsize(i) = mean(volume(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1)));
        zdens(i) = numel([currcent(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1),2)]);  
        ztints(i) = sum(tints(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1)));      
    end
end
normset = hw1; usebins = binedges;
geph_ztints = ztints./(normset'*(dsvoxel*10^-3)^3);
geph_zdens = zdens./(normset'*(dsvoxel*10^-3)^3);

subplot(1,3,1)
bar(usebins*0.001',zsize*0.0158*0.0158*0.07,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
subplot(1,3,2)
bar(usebins*0.001',geph_ztints,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
subplot(1,3,3)
bar(usebins*0.001',geph_zdens,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
savefig([base_folder 'paired_gephyrin_laminar_distributions_withedges.fig'])
%%
figure;
currcent = centPs(~pairedp_idx,:);
currstat = sizeshape_matPs(~pairedp_idx,:);
volume = currstat(:,19);
tints = currstat(:,20);
clear zsize  ztints zdens numbins hy2

for i=1:(numel(binedges))
    if i==(numel(binedges))
        zsize(i) = mean(volume(currcent(:,2)== binedges(i)));
        zdens(i) = numel([currcent(currcent(:,2)== binedges(i),2)]);
        ztints(i) = sum(tints(currcent(:,2)== binedges(i)));
    else
        zsize(i) = mean(volume(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1)));
        zdens(i) = numel([currcent(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1),2)]);  
        ztints(i) = sum(tints(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1)));      
    end
end
normset = hw1; usebins = binedges;
pre_unpaired_ztints = ztints./(normset'*(dsvoxel*10^-3)^3);
pre_unpaired_zdens = zdens./(normset'*(dsvoxel*10^-3)^3);

subplot(1,3,1)
bar(usebins*0.001',zsize*0.0158*0.0158*0.07,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
subplot(1,3,2)
bar(usebins*0.001',pre_unpaired_ztints,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
subplot(1,3,3)
bar(usebins*0.001',pre_unpaired_zdens,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
savefig([base_folder 'unpaired_presynaptic_laminar_distributions_withedges.fig'])

%
figure;
currcent = centPs(pairedp_idx,:);
currstat = sizeshape_matPs(pairedp_idx,:);
volume = currstat(:,19);
tints = currstat(:,20);
clear zsize  ztints zdens numbins hy2

for i=1:(numel(binedges))
    if i==(numel(binedges))
        zsize(i) = mean(volume(currcent(:,2)== binedges(i)));
        zdens(i) = numel([currcent(currcent(:,2)== binedges(i),2)]);
        ztints(i) = sum(tints(currcent(:,2)== binedges(i)));
    else
        zsize(i) = mean(volume(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1)));
        zdens(i) = numel([currcent(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1),2)]);  
        ztints(i) = sum(tints(currcent(:,2)>binedges(i) & ...
            currcent(:,2)<binedges(i+1)));      
    end
end
normset = hw1; usebins = binedges;
pre_paired_ztints = ztints./(normset'*(dsvoxel*10^-3)^3);
pre_paired_zdens = zdens./(normset'*(dsvoxel*10^-3)^3);

subplot(1,3,1)
bar(usebins*0.001',zsize*0.0158*0.0158*0.07,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
subplot(1,3,2)
bar(usebins*0.001',pre_paired_ztints,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
subplot(1,3,3)
bar(usebins*0.001',pre_paired_zdens,'FaceColor',[0 0 1],'BarWidth',1)
set(gca,'XLim',[0 46])
savefig([base_folder 'paired_presynaptic_laminar_distributions_withedges.fig'])

%
i=22;
figure;
bar(usebins(i:86)*0.001',(zscore(pre_paired_ztints(i:86))-...
    zscore(pre_unpaired_ztints(i:86))))
hold on
plot(0:46,zeros(37,1),'k')
set(gca,'XLim',[0 46])
savefig([base_folder 'paired_presynaptic_index_with_edges.fig'])

