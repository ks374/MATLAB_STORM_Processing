clear all

%first separate synaptic groups,for this pass, 
%first  do it quickly with manual separation of synaptic and non_synaptic
%%
%base_folder = 'Y:/backup/ysigal/GLD3_bistratified_RGC/elastic_thresh/';
base_folder = '/n/contefs1/backup/ysigal/GLD3_bistratified_RGC/elastic_thresh/';
base_folder1 = '/n/contefs1/backup/ysigal/GLD3_bistratified_RGC/elastic_align/';
%
load([base_folder 'sizeshapematP_and_cent_water10.mat']);
%
centP = centGw;
%
sizeshape_matP = sizeshape_mat;
voxel = [15.8, 15.8, 70];
for i=1:numel(centP(:,1))
   centP(i,:) = centP(i,:).*voxel ;
end
%
load([base_folder1 'rotated_ds_data.mat']);
%%
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
kk = regionbounds(2,1)-regionbounds(1,1)-1;
lamsize = max(CW(:,1))-min(CW(:,1));
samp_freq = kk;
binedges = 0:lamsize/samp_freq:lamsize;
hw1 = histc(CW(:,1),binedges);%
figure;
bar(binedges*0.001, hw1*((dsvoxel*10^-3)^3))
%
ccent = centP;
for jj = 1:numel(ccent(:,1))   
    ccent(jj,2) = ccent(jj,2)/(dsvoxel);
    ccent(jj,1) = ccent(jj,1)/(dsvoxel);
    ccent(jj,3) = ccent(jj,3)/(dsvoxel);
    [ccent(jj,2), ccent(jj,1), ccent(jj,3)] = ...
        tformfwd(tform2, [ccent(jj,2), ccent(jj,1), ccent(jj,3)]);
    for k=1:3
         ccent(jj,k) = ccent(jj,k)*dsvoxel; 
    end
end
ccent(:,1)=ccent(:,1)-regionbounds(1,2)*dsvoxel;
ccent(:,2)=ccent(:,2)-regionbounds(1,1)*dsvoxel;
ccent(:,3)=ccent(:,3)-regionbounds(1,3)*dsvoxel;
rcentP = ccent;
%
figure;
hist(rcentP(:,2),100)
%
centPwa = centGw(rcentP(:,2)>0,:);
centPa = centP(rcentP(:,2)>0,:);
rcentPa = rcentP(rcentP(:,2)>0,:);
sizeshape_matPa = sizeshape_matP(rcentP(:,2)>0,:);
%
%load([base_folder 'synapse_region_use_edges.mat'])
%for i=1:numel(centPwa(:,1))
%centPwa(i,:) = centPwa(i,:)./[10 10 1];
%end
%in = ~inpolyhedron(r2FV,centPwa);
%numel(find(in))
%
%centPb = centPa(in,:);
%rcentPb = rcentPa(in,:);
%sizeshape_matPb = sizeshape_matPa(in,:);
%
%save([base_folder 'synapse_region_inside_edges_Pw10.mat'],'in')
%%

sizeshapuse = sizeshape_matPa;
val1w = (sizeshapuse(:,14)+1)./(sizeshapuse(:,6)+1);
val2w = log10((sizeshapuse(:,14)+1)*0.0158*0.0158*0.07);
numel(find(val1w<0.99))
%
Xn=80; Yn=70;
Xrange=[0.3 1.0]; Yrange=[-3.5 -0.5];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ; 

figure;
H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
close
%%
cutoffg =700;
H2 = H;
H2(H2>cutoffg)=cutoffg;
%
figure;
pcolor(X,Y,H2)
%%
savefig([base_folder 'Psizeshape_heatmap.fig'])

%savefig([base_folder 'Psizeshape_heatmap_includingedges.fig'])
%%
threshx = zeros(numel(X),1);
threshx1 = zeros(numel(X),1);
%threshx2 = zeros(numel(X),1);
clear all_cvals
%parpool(32)
parfor i=1:numel(X)

test2 = H(:,i);
[xData, yData] = prepareCurveData( Y, test2);

% Set up fittype and options.
ft = fittype( 'gauss2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -4 0.3 0 -2.2 0.3];
opts.StartPoint = [300 -3 0.3 200 -2 0.4];
opts.Upper = [Inf -2 1 cutoffg 0 1];

%ft = fittype( 'gauss3' );
%opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Display = 'Off';
%opts.Lower = [0 -7.3 0.3 0 -8 0.3 0 -6 0.2];
%opts.StartPoint = [100 -7 0.2 100 -7 0.2 100 -5 0.2];
%opts.Upper = [Inf -6.8 1 Inf -6.8 2 cutoff 0 2];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult, xData, yData );
xData2 = -4:0.01:0;
coeffvals = coeffvalues(fitresult);
all_cvals(i,:) = coeffvals;
y1 = coeffvals(1)*exp(-((xData2-coeffvals(2))/coeffvals(3)).^2);
y2 =coeffvals(4)*exp(-((xData2-coeffvals(5))/coeffvals(6)).^2);
threshx(i) = xData2(find(cumsum(y1/sum(y1))>0.75,1,'first'));

end
%
%
figure;
plot(X,threshx,'r*')


alpha(0.3)
hold on
pcolor(X,Y,H2)
%%
all_cvals
X2 = X;
threshx2 = threshx;
X2(threshx2==0)=[];
threshx2(threshx2==0)=[];

%X2(all_cvals(:,5)<-5.95)=[];
%threshx2(all_cvals(:,5)<-5.95)=[];
%
X2(all_cvals(:,1)<100 | all_cvals(:,4)<100 | all_cvals(:,5)>-2)=[];
threshx2(all_cvals(:,1)<100 | all_cvals(:,4)<100 | all_cvals(:,5)>-2)=[];
%
[xData, yData] = prepareCurveData( X2, threshx2 );

% Set up fittype and options.
ft = fittype( 'poly3' );
clear fitresultg coeffvals2
% Fit model to data.
[fitresultg, gof] = fit( xData, yData, ft );
coeffvals2 = coeffvalues(fitresultg);
xfit = 0.3:.01:1;
%yfit = coeffvals2(1)*xfit.^4 + coeffvals2(2)*xfit.^3 + ...
%    coeffvals2(3)*xfit.^2 + coeffvals2(4)*xfit + coeffvals2(5);
yfit = coeffvals2(1)*xfit.^3 + coeffvals2(2)*xfit.^2 + ...
    coeffvals2(3)*xfit + coeffvals2(4);
%yfit = coeffvals2(1)*xfit.^2 + ...
%    coeffvals2(2)*xfit + coeffvals2(3);

% Plot fit with data.
figure;
plot( xData, yData,'r.','MarkerSize',20 );
hold on

plot(xfit,yfit,'r','LineWidth',5)
%
alpha(0.3)
hold on
pcolor(X,Y,H2)
%
%figure;
%plot(gofall)
%savefig([base_folder 'Pfit_sizeshape_heatmap.fig'])
%
yvals = coeffvals2(1)*val1w.^3 + coeffvals2(2)*val1w.^2 + ...
    coeffvals2(3)*val1w + coeffvals2(4);
%yvals = coeffvals2(1)*val1.^2 + coeffvals2(2)*val1 + ...
%    coeffvals2(3);
idxPselect = val2w>=yvals;
numel(find(idxPselect))

%idxPselect2 = val2w>=yvals;
numel(find(idxPselect2))
%%
load([base_folder 'statsPwater10.mat']);
%%
statsPa = statsGwater(rcentP(:,2)>0,:);
%%

centPa2s = centPa(idxPselect2,:);
rcentPa2s = rcentPa(idxPselect2,:);
sizeshape_matPa2s = sizeshape_matPa(idxPselect2,:);
centPa2ns = centPa(~idxPselect2,:);
rcentPa2ns = rcentPa(~idxPselect2,:);
sizeshape_matPa2ns = sizeshape_matPa(~idxPselect2,:);
%
statsPa2s = statsPa(idxPselect2,:);
statsPa2ns = statsPa(~idxPselect2,:);
centPa2ns = centPa2ns(sizeshape_matPa2ns(:,19)>50,:);
rcentPa2ns = rcentPa2ns(sizeshape_matPa2ns(:,19)>50,:);
statsPa2ns = statsPa2ns(sizeshape_matPa2ns(:,19)>50,:);
sizeshape_matPa2ns = sizeshape_matPa2ns(sizeshape_matPa2ns(:,19)>50,:);
%
centPa1s = centPa(idxPselect,:);
rcentPa1s = rcentPa(idxPselect,:);
sizeshape_matPa1s = sizeshape_matPa(idxPselect,:);
centPa1ns = centPa(~idxPselect,:);
rcentPa1ns = rcentPa(~idxPselect,:);
sizeshape_matPa1ns = sizeshape_matPa(~idxPselect,:);
%
statsPa1s = statsPa(idxPselect,:);
statsPa1ns = statsPa(~idxPselect,:);
centPa1ns = centPa1ns(sizeshape_matPa1ns(:,19)>50,:);
rcentPa1ns = rcentPa1ns(sizeshape_matPa1ns(:,19)>50,:);
statsPa1ns = statsPa1ns(sizeshape_matPa1ns(:,19)>50,:);
sizeshape_matPa1ns = sizeshape_matPa1ns(sizeshape_matPa1ns(:,19)>50,:);
%
save([base_folder 'statslistP2sw10_lo_thresh.mat'],'centPa1s','rcentPa1s','sizeshape_matPa1s','-v7.3')
save([base_folder 'statslistP2nsw10_lo_thresh.mat'],'centPa1ns','rcentPa1ns','sizeshape_matPa1ns','-v7.3')
%
save([base_folder 'statslistP2sw10.mat'],'centPa2s','rcentPa2s','sizeshape_matPa2s','-v7.3')
save([base_folder 'statslistP2nsw10.mat'],'centPa2ns','rcentPa2ns','sizeshape_matPa2ns','-v7.3')
%
save([base_folder 'statsP2sw10_lo_thresh.mat'],'statsPa1s','centPa1s','rcentPa1s','sizeshape_matPa1s','-v7.3')
save([base_folder 'statsP2nsw10_lo_thresh.mat'],'statsPa1ns','centPa1ns','rcentPa1ns','sizeshape_matPa1ns','-v7.3')
save([base_folder 'statsP2sw10.mat'],'statsPa2s','centPa2s','rcentPa2s','sizeshape_matPa2s','-v7.3')
save([base_folder 'statsP2nsw10.mat'],'statsPa2ns','centPa2ns','rcentPa2ns','sizeshape_matPa2ns','-v7.3')

