%%
outpath = 'Z:\Jackie\neuroD6_retina_round2_section10_3\9_4_2019_Processing\';
load([outpath 'sizeshapematR_and_cent_water10_area_cutoff.mat']);
load([outpath 'rotated_ds_data.mat']);
dsvoxel  =158/0.3;
%%
centG = centGw;
sizeshape_matG = sizeshape_mat;
voxel = [15.8, 15.8, 70];
for i=1:numel(centG(:,1))
   centG(i,:) = centG(i,:).*voxel ;
end
%%
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
%
ccent = centG;
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
rcentG = ccent;
%
centGwa = centGw(rcentG(:,2)>0,:);
centGa = centG(rcentG(:,2)>0,:);
rcentGa = rcentG(rcentG(:,2)>0,:);
sizeshape_matGa = sizeshape_matG(rcentG(:,2)>0,:);

%%
i=6;  % i = range(8)
sizeshapuse = sizeshape_matGa;
val1w = (sizeshapuse(:,10+i)+1)./(sizeshapuse(:,2+i)+1);
val2w = log10((sizeshapuse(:,10+i)+1)*0.0158*0.0158*0.07);
numel(find(val1w<0.99))
%
Xn=80; Yn=60;
Xrange=[0.4 1.0]; Yrange=[-3.5 -0.5];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ; 

figure;
H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
close

cutoffg =60;
H2 = H;
H2(H2>cutoffg)=cutoffg;

figure;
pcolor(X,Y,H2)

savefig([outpath 'Rsizeshape_heatmap.fig'])
%%
threshx = zeros(numel(X),1);
threshx1 = zeros(numel(X),1);
%threshx2 = zeros(numel(X),1);
clear all_cvals

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
threshx(i) = xData2(find(cumsum(y1/sum(y1))>0.5,1,'first'));

end
%
%%
figure;
pcolor(X,Y,H2)
alpha(0.8)
hold on
plot(X,threshx,'r*')
%%
X2 = X;
threshx2 = threshx;
X2(threshx2==0)=[];
threshx2(threshx2==0)=[];

%X2(all_cvals(:,5)<-5.95)=[];
%threshx2(all_cvals(:,5)<-5.95)=[];
%
% X2(all_cvals(:,1)<40 | all_cvals(:,4)<60 | all_cvals(:,6)<0.1)=[];
% threshx2(all_cvals(:,1)<40 | all_cvals(:,4)<60 | all_cvals(:,6)<0.1)=[];
%
[xData, yData] = prepareCurveData( X2, threshx2 );

% Set up fittype and options.
ft = fittype( 'poly2' );
clear fitresultg coeffvals2
% Fit model to data.
[fitresultg, gof] = fit( xData, yData, ft );
coeffvals2 = coeffvalues(fitresultg);
xfit = 0.40:.01:1;
%yfit = coeffvals2(1)*xfit.^4 + coeffvals2(2)*xfit.^3 + ...
%    coeffvals2(3)*xfit.^2 + coeffvals2(4)*xfit + coeffvals2(5);
%yfit = coeffvals2(1)*xfit.^3 + coeffvals2(2)*xfit.^2 + ...
%    coeffvals2(3)*xfit + coeffvals2(4);
yfit = coeffvals2(1)*xfit.^2 + ...
    coeffvals2(2)*xfit + coeffvals2(3);

% Plot fit with data.
figure;
pcolor(X,Y,H2)
hold on
alpha(0.5)
plot( xData, yData,'b.','MarkerSize',20 );
hold on

plot(xfit,yfit,'r','LineWidth',5)

%%
% figure;
% plot(gofall)
savefig([outpath 'Rfit_sizeshape_heatmap.fig'])
%%
%yvals = coeffvals2(1)*val1w.^3 + coeffvals2(2)*val1w.^2 + ...
%    coeffvals2(3)*val1w + coeffvals2(4);
yvals = coeffvals2(1)*val1w.^2 + coeffvals2(2)*val1w + ...
    coeffvals2(3);
%idxGselect = val2w>=yvals;
%numel(find(idxGselect))

idxGselect2 = val2w>=yvals;
numel(find(idxGselect2))
%%
load([outpath 'statsRwater10_area_cutoff.mat']);
%
statsGa = statsGwater(rcentG(:,2)>0,:);
%statsGb = statsGa(in,:);
%
%centGb2s = centGb(idxGselect,:);
%rcentGb2s = rcentGb(idxGselect,:);
%sizeshape_matGb2s = sizeshape_matGb(idxGselect,:);
%centGb2ns = centGb(~idxGselect,:);
%rcentGb2ns = rcentGb(~idxGselect,:);
%sizeshape_matGb2ns = sizeshape_matGb(~idxGselect,:);
%statsGb2s = statsGb(idxGselect,:);
%statsGb2ns = statsGb(~idxGselect,:);

%centGb2ns = centGb2ns(sizeshape_matGb2ns(:,19)>50,:);
%rcentGb2ns = rcentGb2ns(sizeshape_matGb2ns(:,19)>50,:);
%statsGb2ns = statsGb2ns(sizeshape_matGb2ns(:,19)>50,:);
%sizeshape_matGb2ns = sizeshape_matGb2ns(sizeshape_matGb2ns(:,19)>50,:);
%
centGa2s = centGa(idxGselect2,:);
rcentGa2s = rcentGa(idxGselect2,:);
sizeshape_matGa2s = sizeshape_matGa(idxGselect2,:);
centGa2ns = centGa(~idxGselect2,:);
rcentGa2ns = rcentGa(~idxGselect2,:);
sizeshape_matGa2ns = sizeshape_matGa(~idxGselect2,:);
%
statsGa2s = statsGa(idxGselect2,:);
statsGa2ns = statsGa(~idxGselect2,:);
%
centGa2ns = centGa2ns(sizeshape_matGa2ns(:,19)>50,:);
rcentGa2ns = rcentGa2ns(sizeshape_matGa2ns(:,19)>50,:);
%statsGa2ns = statsGa2ns(sizeshape_matGa2ns(:,19)>50,:);
sizeshape_matGa2ns = sizeshape_matGa2ns(sizeshape_matGa2ns(:,19)>50,:);
%%
%save([base_folder 'statslistG2sw10.mat'],'centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statslistG2nsw10.mat'],'centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
save([outpath 'statslistR2sw10.mat'],'centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statslistR2nsw10.mat'],'centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')

%
save([outpath 'statsR2sw10.mat'],'statsGa2s','centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statsR2nsw10.mat'],'statsGa2ns','centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')
%save([base_folder 'statsG2sw10.mat'],'statsGb2s','centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statsG2nsw10.mat'],'statsGb2ns','centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
