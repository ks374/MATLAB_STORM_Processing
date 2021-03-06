%%
outpath = 'Z:\Chenghang\chenghaz_008_Sample_Epi_Ipsi_A\analysis\elastic_align\Result\';
load([outpath 'sizeshapematG_and_cent_water10_area_cutoff.mat']);
dsvoxel  =155/0.3;
%
centG = centGw;
sizeshape_matG = sizeshape_mat;
voxel = [15.5, 15.5, 70];
for i=1:numel(centG(:,1))
   centG(i,:) = centG(i,:).*voxel ;
end
rcentG = centG;
%%
i=4;  % i = range(8)
sizeshapuse = sizeshape_mat;
val1w = (sizeshapuse(:,10+i)+1)./(sizeshapuse(:,2+i)+1);
val2w = log10((sizeshapuse(:,10+i)+1)*0.0158*0.0158*0.07);
numel(find(val1w<0.99))

Xn=80; Yn=60;
Xrange=[0.4 1.0]; Yrange=[-3.5 -0.5];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ; 

figure;
H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange);
close

cutoffg = 1;
H2 = H;
H2(H2>cutoffg)=cutoffg;

figure;
pcolor(X,Y,H2)
%%
savefig([outpath 'Gsizeshape_heatmap.fig'])
%%
load([outpath 'statsGwater10_area_cutoff.mat']);
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
%%
% threshx = zeros(numel(X),1);
% threshx1 = zeros(numel(X),1);
% %threshx2 = zeros(numel(X),1);
% clear all_cvals
% 
% parfor i=1:numel(X)
% 
% test2 = H(:,i);
% [xData, yData] = prepareCurveData( Y, test2);
% 
% % Set up fittype and options.
% ft = fittype( 'gauss2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [0 -4 0.3 0 -2.2 0.3];
% opts.StartPoint = [300 -3 0.3 200 -2 0.4];
% opts.Upper = [Inf -2 1 cutoffg 0 1];
% 
% %ft = fittype( 'gauss3' );
% %opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% %opts.Display = 'Off';
% %opts.Lower = [0 -7.3 0.3 0 -8 0.3 0 -6 0.2];
% %opts.StartPoint = [100 -7 0.2 100 -7 0.2 100 -5 0.2];
% %opts.Upper = [Inf -6.8 1 Inf -6.8 2 cutoff 0 2];
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% %figure( 'Name', 'untitled fit 1' );
% %h = plot( fitresult, xData, yData );
% xData2 = -4:0.01:0;
% coeffvals = coeffvalues(fitresult);
% all_cvals(i,:) = coeffvals;
% y1 = coeffvals(1)*exp(-((xData2-coeffvals(2))/coeffvals(3)).^2);
% y2 =coeffvals(4)*exp(-((xData2-coeffvals(5))/coeffvals(6)).^2);
% threshx(i) = xData2(find(cumsum(y1/sum(y1))>0.7,1,'first'));
% 
% end
% %
% 
% figure;
% pcolor(X,Y,H2)
% alpha(0.8)
% hold on
% plot(X,threshx,'r*')
% %%
% X2 = X;
% threshx2 = threshx;
% X2(threshx2==0)=[];
% threshx2(threshx2==0)=[];
% 
% %X2(all_cvals(:,5)<-5.95)=[];
% %threshx2(all_cvals(:,5)<-5.95)=[];
% %
% % X2(all_cvals(:,1)<40 | all_cvals(:,4)<60 | all_cvals(:,6)<0.1)=[];
% % threshx2(all_cvals(:,1)<40 | all_cvals(:,4)<60 | all_cvals(:,6)<0.1)=[];
% %
% [xData, yData] = prepareCurveData( X2, threshx2 );
% 
% % Set up fittype and options.
% ft = fittype( 'smoothingspline');
% clear fitresultg coeffvals2
% % Fit model to data.
% [fitresultg, gof] = fit( xData, yData, ft , 'SmoothingParam',0.9999999999 );
% 
% xfit = 0.30:.01:1;
% %yfit = coeffvals2(1)*xfit.^4 + coeffvals2(2)*xfit.^3 + ...
% %    coeffvals2(3)*xfit.^2 + coeffvals2(4)*xfit + coeffvals2(5);
% %yfit = coeffvals2(1)*xfit.^3 + coeffvals2(2)*xfit.^2 + ...
% %    coeffvals2(3)*xfit + coeffvals2(4);
% yfit = fitresultg(xfit);
% 
% % Plot fit with data.
% figure;
% pcolor(X,Y,H2)
% hold on
% alpha(0.5)
% plot( xData, yData,'b.','MarkerSize',20 );
% hold on
% 
% plot(xfit,yfit,'r','LineWidth',5)
%%
% X2 = X;
% threshx2 = threshx;
% X2(threshx2==0)=[];
% threshx2(threshx2==0)=[];
% 
% %X2(all_cvals(:,5)<-5.95)=[];
% %threshx2(all_cvals(:,5)<-5.95)=[];
% %
% % X2(all_cvals(:,1)<40 | all_cvals(:,4)<60 | all_cvals(:,6)<0.1)=[];
% % threshx2(all_cvals(:,1)<40 | all_cvals(:,4)<60 | all_cvals(:,6)<0.1)=[];
% %
% [xData, yData] = prepareCurveData( X2, threshx2 );
% 
% % Set up fittype and options.
% ft = fittype( 'poly2' );
% clear fitresultg coeffvals2
% % Fit model to data.
% [fitresultg, gof] = fit( xData, yData, ft );
% coeffvals2 = coeffvalues(fitresultg);
% xfit = 0.40:.01:1;
% %yfit = coeffvals2(1)*xfit.^4 + coeffvals2(2)*xfit.^3 + ...
% %    coeffvals2(3)*xfit.^2 + coeffvals2(4)*xfit + coeffvals2(5);
% %yfit = coeffvals2(1)*xfit.^3 + coeffvals2(2)*xfit.^2 + ...
% %    coeffvals2(3)*xfit + coeffvals2(4);
% yfit = coeffvals2(1)*xfit.^2 + ...
%     coeffvals2(2)*xfit + coeffvals2(3);
% 
% % Plot fit with data.
% figure;
% pcolor(X,Y,H2)
% hold on
% alpha(0.5)
% plot( xData, yData,'b.','MarkerSize',20 );
% hold on
% 
% plot(xfit,yfit,'r','LineWidth',5)

%%
%figure;
%plot(gofall)
% savefig([outpath 'Gfit_sizeshape_heatmap.fig'])
%%
% %yvals = coeffvals2(1)*val1w.^3 + coeffvals2(2)*val1w.^2 + ...
% %    coeffvals2(3)*val1w + coeffvals2(4);
% % yvals = coeffvals2(1)*val1w.^2 + coeffvals2(2)*val1w + ...
% %     coeffvals2(3);
% yvals = fitresultg(val1w);
% %idxGselect = val2w>=yvals;
% %numel(find(idxGselect))
% 
% idxGselect2 = val2w>=yvals;
% numel(find(idxGselect2))
%%
% load([outpath 'statsGwater10_area_cutoff.mat']);
% %
% statsGa = statsGwater(rcentG(:,2)>0,:);
% centGa = centG;
% rcentGa = rcentG;
% sizeshape_matGa = sizeshape_matG;
% %statsGb = statsGa(in,:);
% %
% %centGb2s = centGb(idxGselect,:);
% %rcentGb2s = rcentGb(idxGselect,:);
% %sizeshape_matGb2s = sizeshape_matGb(idxGselect,:);
% %centGb2ns = centGb(~idxGselect,:);
% %rcentGb2ns = rcentGb(~idxGselect,:);
% %sizeshape_matGb2ns = sizeshape_matGb(~idxGselect,:);
% %statsGb2s = statsGb(idxGselect,:);
% %statsGb2ns = statsGb(~idxGselect,:);
% 
% %centGb2ns = centGb2ns(sizeshape_matGb2ns(:,19)>50,:);
% %rcentGb2ns = rcentGb2ns(sizeshape_matGb2ns(:,19)>50,:);
% %statsGb2ns = statsGb2ns(sizeshape_matGb2ns(:,19)>50,:);
% %sizeshape_matGb2ns = sizeshape_matGb2ns(sizeshape_matGb2ns(:,19)>50,:);
% %
% centGa2s = centGa(idxGselect2,:);
% rcentGa2s = rcentGa(idxGselect2,:);
% sizeshape_matGa2s = sizeshape_matGa(idxGselect2,:);
% centGa2ns = centGa(~idxGselect2,:);
% rcentGa2ns = rcentGa(~idxGselect2,:);
% sizeshape_matGa2ns = sizeshape_matGa(~idxGselect2,:);
% %
% statsGa2s = statsGa(idxGselect2,:);
% statsGa2ns = statsGa(~idxGselect2,:);
%%
Real_non_select_area_thre = 10;
centGa2ns = centGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
rcentGa2ns = rcentGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
statsGa2ns = statsGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);
sizeshape_matGa2ns = sizeshape_matGa2ns(sizeshape_matGa2ns(:,19)>Real_non_select_area_thre,:);

%save([base_folder 'statslistG2sw10.mat'],'centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statslistG2nsw10.mat'],'centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
save([outpath 'statslistG2sw10.mat'],'centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statslistG2nsw10.mat'],'centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')

%
save([outpath 'statsG2sw10.mat'],'statsGa2s','centGa2s','rcentGa2s','sizeshape_matGa2s','-v7.3')
save([outpath 'statsG2nsw10.mat'],'statsGa2ns','centGa2ns','rcentGa2ns','sizeshape_matGa2ns','-v7.3')
%save([base_folder 'statsG2sw10.mat'],'statsGb2s','centGb2s','rcentGb2s','sizeshape_matGb2s','-v7.3')
%save([base_folder 'statsG2nsw10.mat'],'statsGb2ns','centGb2ns','rcentGb2ns','sizeshape_matGb2ns','-v7.3')
