%normalize and threshold conventional images after selecting IPL
clear all
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
%%
%path definition: 
conv_path = [exp_folder  'conv_merged_ds/'];
stormfiles = dir([conv_path '*.tif']);
num_images = numel(stormfiles);e
load([exp_folder 'rotated_ds_data.mat']);
%%
num_images2 = size(BG4,3);
sec_use = 1:num_images;
sec_use(badsec)=[];
%%
%input: num_images2, conv_path, stormfiles, 

hy1c = zeros(num_images2,255);
hy2c = zeros(num_images2,255);
hy3c = zeros(num_images2,255);

parfor j=1:num_images2
    disp(strcat('processing image #',sprintf('%03d',j))); 
    A = (imread([conv_path stormfiles(sec_use(j),1).name]));
    A1 = A(:,:,1);
    A2 = A(:,:,2);
%    A3 = A(:,:,3);
    B1 = imresize(BP4(:,:,j),[size(A1,1) size(A1,2)]);
    B1 = logical(B1);
    A1 = A1.*uint8(B1);
    A2 = A2.*uint8(B1);
    %A3 = A3.*uint8(B1);
    %figure;
    %imshow(A1)
    A1a = A1(find(A1));
    A2a = A2(find(A2));
    %A3a = A3(find(A3));
    h_all(j) = numel(find(B1));
    [hy hx] = hist(A1a,1:1:255);
    hy1c(j,:) = hy;
    [hy hx] = hist(A2a,1:1:255);
    hy2c(j,:) = hy;
    %[hy hx] = hist(A3a,1:1:255);
    %hy3c(j,:) = hy;
    %hy3(j,:) = hy./h_all(j);
end
%
x_hist = 1:255;
x_sec = 1:num_images2;
%
h_all2 = h_all./mean(h_all);
for i=1:num_images2
    hy1cb(i,:) = hy1c(i,:)./h_all2(i);
    hy2cb(i,:) = hy2c(i,:)./h_all2(i);
%    hy3cb(i,:) = hy3c(i,:)./sum(hy3c(i,:));
%
%    hy1_b(i,:) = hy1(i,:)./sum(hy1(i,:));
%    hy2_b(i,:) = hy2(i,:)./sum(hy2(i,:));
%    hy3_b(i,:) = hy3(i,:)./sum(hy3(i,:));
end
%
chan = hy1cb;
for idx = 1:255
    disp(idx)
    zthresh = 3;
    currparam = chan(:,idx);
currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
%[xData, yData, weights] = prepareCurveData( currfitx', currfity, 5*(717-currfitx)' );
[xData, yData] = prepareCurveData( currfitx', currfity );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
%opts.SmoothingParam = 0.9999999646249706;
%opts.Weights = weights;

opts.SmoothingParam = 0.7;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult,  xData',  yData' );
varuse1(idx,:) = feval(fitresult,1:num_images);
end
%
chan = hy2cb;
for idx = 1:255
    disp(idx)
    zthresh = 3;
    currparam = chan(:,idx);
currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
%[xData, yData, weights] = prepareCurveData( currfitx', currfity, 5*(717-currfitx)' );
[xData, yData] = prepareCurveData( currfitx', currfity );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
%opts.SmoothingParam = 0.9999999646249706;
%opts.Weights = weights;

opts.SmoothingParam = 0.7;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult,  xData',  yData' );
varuse2(idx,:) = feval(fitresult,1:num_images);
end
%%
figure;
plot(varuse1(:,100))
%%
path4 = [exp_folder  'elastic_norm/'];
if exist(path4)~=7
    mkdir(path4);
	mkdir([path4 'storm_merged/'])
	mkdir([path4 'storm_merged_ds/'])
	mkdir([path4 'conv_merged/'])
	mkdir([path4 'conv_merged_ds/'])	
    mkdir([path4 'conv_561/'])
	mkdir([path4 'conv_561_ds/'])
end
%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
%%
figure;
plot(hgram2)
hold on 
plot(hy2c(i,:),'r')

    B1a = B2(find(B2));
    %A3a = A3(find(A3));
    [hy hx] = hist(B1a,1:1:255);
    hold on 
    plot(hy,'g')

%%
info = imfinfo([conv_path stormfiles(sec_use(1),1).name]);
    width = info.Width;
    height = info.Height;
parfor i=1:num_images2
    disp(i)
hgram1 = varuse1(:,i).*h_all2(i);
hgram2 = varuse2(:,i).*h_all2(i);
%hgram3 = varuse3(:,i)/sum(varuse3(:,i));
%
%figure;
%plot(hgram1)
A = (imread([conv_path stormfiles(sec_use(i),1).name]));
    A1 = A(:,:,1);
    A2 = A(:,:,2);
    A3 = A(:,:,3);

    C1 = imresize(BP4(:,:,i),[size(A1,1) size(A1,2)]);
    C1 = logical(C1);
    A1 = A1.*uint8(C1);
    A2 = A2.*uint8(C1);
    A3 = A3.*uint8(C1);
    %rescaleP(i) = avgP(i)/meanp; 
    %B1 = uint8(double(A1)./rescaleP(i)); 
    %figure;
    %imshow(A1)
    %rescaleG(i) = avgG(i)/meang; 
    %B2 = uint8(double(A2)./rescaleG(i));
     
    
hgram1a = cat(1,(round(width*height-sum(hgram1))),hgram1);
hgram2a = cat(1,(round(width*height-sum(hgram2))),hgram2);
%hgram3a = cat(1,numel(find(A3==0)),hgram3*numel(find(A3)));
%
%figure;
%plot(hgram1a)
B1 = histeq(A1,hgram1a);
%figure;
%imshow(B1)
B2 = histeq(A2,hgram2a);
%B3 = histeq(A3,hgram3a);

%B1=0;
B1(A1<1)=0;
B2(A2<1)=0;
%B3(A3<1)=0;
C = cat(3,B1,B2,A3);
%    A1a = B1(find(B1));
%    h_all(j) = numel(find(B1));
%    [hy hx] = hist(A1a,1:1:255);
%    hy2 = hy/h_all(j);
imwrite(C, [path4 'conv_merged_ds/' ...
        sprintf('%03d',i) '.png']);

%
%figure;
%plot(hy1c(i,:))
%hold on
%plot(hy,'.')
%hold on 
%plot(varuse1(:,i)*numel(find(A1)),'r')
end
%%
path5 = [path4 'conv_merged_ds/'];
files2 = dir([path5 '*.png']);
info = imfinfo([path5 files2(1,1).name]);
    width = info.Width;
    height = info.Height;
A1 = zeros(height,width, num_images2,'uint8'); A2= A1;
parfor i=1:num_images2
A = (imread([path4 'conv_merged_ds/' ...
        sprintf('%03d',i) '.png']));
    A1(:,:,i) = A(:,:,1);
    A2(:,:,i) = A(:,:,2);
    %A3 = A(:,:,3);
end
%%
B1 = logical(A1);
C1 = squeeze(sum(B1,1))+1;
temp1 = squeeze(mean(A1,1))./C1;
B2 = logical(A2);
C2 = squeeze(sum(B2,1))+1;
temp2 = squeeze(mean(A2,1))./C2;

figure;
imshow(uint8(temp1*713*6))
figure;
imshow(uint8(temp2*713*20))
    %%
%now determine threshold based on ds conv image stack
path5 = [path4 'conv_merged_ds/'];
clear B*
files2 = dir([path5 '*.png']);

    info = imfinfo([path5 files2(1,1).name]);
    width = info.Width;
    height = info.Height;
B1 = zeros(1000000000,1); B2 = B1; B3 = B1;
idx1 = 1; idx2 = idx1; idx3=idx1;
for i=1:num_images2
    disp(i)
    A = (imread([path5 files2(i,1).name]));
    x = A(:,:,1);
    B1(idx1:idx1 + numel(find(A(:,:,1)))-1) =  x(find(x));
    x = A(:,:,2);
    B2(idx2:idx2 + numel(find(A(:,:,2)))-1) =  x(find(x));
    x = A(:,:,3);
    B3(idx3:idx3 + numel(find(A(:,:,3)))-1) =  x(find(x));

    idx1 = idx1 + numel(find(A(:,:,1)));
    idx2 = idx2 + numel(find(A(:,:,2)));
    idx3 = idx3 + numel(find(A(:,:,3)));
end
B1 = B1(1:idx1-1); B2 = B2(1:idx2-1); B3 = B3(1:idx3-1);
thresh1 = multithresh(B1,2);
thresh2 = multithresh(B2,2);
thresh3 = multithresh(B3,2);

%
path8 = [exp_folder  'elastic_thresh/'];
if exist(path8)~=7
    mkdir(path8);
	mkdir([path8 'storm_merged/'])
	mkdir([path8 'storm_merged_ds/'])
	mkdir([path8 'conv_merged/'])
	mkdir([path8 'conv_merged_ds/'])	
    mkdir([path8 'conv_561/'])
	mkdir([path8 'conv_561_ds/'])
    mkdir([path8 'for_align/'])
	mkdir([path8 'for_align_ds/'])
end
%hy1 = zeros(num_images2,255);
%hy2 = zeros(num_images2,255);
%hy3 = zeros(num_images2,255);
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end
%%
path = [exp_folder 'storm_merged/'];
parfor j=1:num_images2
    disp(j)
     A = (imread([path stormfiles(sec_use(j),1).name]));
     A1 = A(:,:,1);
     A2 = A(:,:,2);
     A3 = A(:,:,3);
 
    A = (imread([path5 files2(j,1).name]));
    A1c = A(:,:,1);
    A2c = A(:,:,2);
    A3c = A(:,:,3);

    A1cbw  = im2bw(A1c,thresh1(1)/256);
    A2cbw  = im2bw(A2c,thresh2(1)/256);
    A3cbw  = im2bw(A3c,thresh3(1)/256);

    A1ct = A1c.*uint8(A1cbw);
    A2ct = A2c.*uint8(A2cbw);
    A3ct = A3c.*uint8(A3cbw);
    
    %A1ct = imadjust((A1ct),[plow phi],[0.0 1]);
    %A2ct = imadjust((A2ct),[glow ghi],[0.0 1]);

    
    A1t = A1.*imresize(uint8(A1cbw),[size(A1,1) size(A1,2)],'nearest');
    A2t = A2.*imresize(uint8(A2cbw),[size(A1,1) size(A1,2)],'nearest');
    A3t = A3.*imresize(uint8(A3cbw),[size(A1,1) size(A1,2)],'nearest');
    
    
    %A1t = imadjust(A1t,stretchlim(A1t,[0 1]),[]);
    %A2t = imadjust(A2t,stretchlim(A2t,[0 1]),[]);
    
     A1a = A1t(find(A1t));
     A2a = A2t(find(A2t));
     A3a = A3t(find(A3t));
     [hy hx] = hist(A1a,1:1:255);
     hy1(j,:) = hy;
     [hy hx] = hist(A2a,1:1:255);
     hy2(j,:) = hy;
     [hy hx] = hist(A3a,1:1:255);
     hy3(j,:) = hy;
C = cat(3,A1ct,A2ct,A3ct);
%figure;
%imshow(C)
imwrite(C, [path8 'conv_merged_ds/' ...
        sprintf('%03d',j) '.png']);
C1 = cat(3,A1t,A2t,A3t);

imwrite(C1, [path8 'storm_merged/' ...
        sprintf('%03d',j) '.png']);
    C1_small = imresize(C1, 0.1);
    imwrite(C1_small, [path8 'storm_merged_ds/' ...
        sprintf('%03d',j) '.png']);
end
%
save([exp_folder 'stormhists.mat'],'hy1','hy2','hy3')
save([exp_folder 'convhists.mat'],'hy1c','hy2c','hy3c','h_all')
%%
load([exp_folder 'stormhists.mat'],'hy1','hy2','hy3')
load([exp_folder 'convhists.mat'],'hy1c','hy2c','hy3c','h_all')
%%
%
x_hist = 1:255;
x_sec = 1:num_images2;
%
for i=1:num_images2
    hy1_s(i,:) = hy1(i,:)./h_all2(i);
    hy2_s(i,:) = hy2(i,:)./h_all2(i);
%    hy3cb(i,:) = hy3c(i,:)./sum(hy3c(i,:));
end
%
chan = hy1_s;
for idx = 1:255
    disp(idx)
    zthresh = 3;
    currparam = chan(:,idx);
currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
%[xData, yData, weights] = prepareCurveData( currfitx', currfity, 5*(717-currfitx)' );
[xData, yData] = prepareCurveData( currfitx', currfity );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
%opts.SmoothingParam = 0.9999999646249706;
%opts.Weights = weights;

opts.SmoothingParam = 0.7;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult,  xData',  yData' );
varuse1s(idx,:) = feval(fitresult,1:num_images);
end
%%
chan = hy2_s;
for idx = 1:255
    disp(idx)
    zthresh = 3;
    currparam = chan(:,idx);
currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
%[xData, yData, weights] = prepareCurveData( currfitx', currfity, 5*(717-currfitx)' );
[xData, yData] = prepareCurveData( currfitx', currfity );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
%opts.SmoothingParam = 0.9999999646249706;
%opts.Weights = weights;

opts.SmoothingParam = 0.1;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult,  xData',  yData' );
varuse2s(idx,:) = feval(fitresult,1:num_images);
end
%%
chan = hy3_s;
for idx = 1:255
    disp(idx)
    zthresh = 3;
    currparam = chan(:,idx);
currfitx = x_sec(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
currfity = currparam(~(zscore(currparam) <-zthresh | zscore(currparam) >zthresh));
%[xData, yData, weights] = prepareCurveData( currfitx', currfity, 5*(717-currfitx)' );
[xData, yData] = prepareCurveData( currfitx', currfity );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.Normalize = 'on';
%opts.SmoothingParam = 0.9999999646249706;
%opts.Weights = weights;

opts.SmoothingParam = 0.4;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult,  xData',  yData' );
varuse3(idx,:) = feval(fitresult,1:num_images2);
end
%%
path9 = [base_path 'GLD3_bistratified_RGC/elastic_thresh/'];
if exist(path9)~=7
    mkdir(path9);
	mkdir([path9 'storm_merged/'])
	mkdir([path9 'storm_merged_ds/'])
	mkdir([path9 'conv_merged/'])
	mkdir([path9 'conv_merged_ds/'])	
    mkdir([path9 'conv_561/'])
	mkdir([path9 'conv_561_ds/'])
end
%
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(32)
end

parfor i=1:num_images2
    disp(i)
hgram1 = varuse1s(:,i).*h_all2(i);
hgram2 = varuse2s(:,i).*h_all2(i);
%
A = (imread([path8 'storm_merged/' ...
        sprintf('%03d',i) '.png']));
    A1 = A(:,:,1);
    A2 = A(:,:,2);
    A3 = A(:,:,3);

hgram1a = cat(1,(round(100*width*height-sum(hgram1))),hgram1);
hgram2a = cat(1,(round(100*width*height-sum(hgram2))),hgram2);

B1 = histeq(A1,hgram1a);
%figure;
%imshow(B1)
B2 = histeq(A2,hgram2a);
%B3 = histeq(A3,hgram3a);

%B1=0;
B1(A1<1)=0;
B2(A2<1)=0;
%B3(A3<1)=0;
C = cat(3,B1,B2,A3);
     
%figure;
%imshow(C)
imwrite(C, [path9 'storm_merged/' ...
        sprintf('%03d',i) '.png']);
    C_small = imresize(C, 0.1);
    imwrite(C_small, [path9 'storm_merged_ds/' ...
        sprintf('%03d',i) '.png']);

end
%%
path5 = [path9 'storm_merged_ds/'];
files2 = dir([path5 '*.png']);
info = imfinfo([path5 files2(1,1).name]);
    width = info.Width;
    height = info.Height;
A1 = zeros(height,width, num_images2,'uint8'); A2= A1;
parfor i=1:num_images2
A = (imread([path9 'storm_merged_ds/' ...
        sprintf('%03d',i) '.png']));
    A1(:,:,i) = A(:,:,1);
    A2(:,:,i) = A(:,:,2);
    %A3 = A(:,:,3);
end
%
B1 = logical(A1);
C1 = squeeze(sum(B1,1))+1;
temp1 = squeeze(mean(A1,1))./C1;
temp1a = squeeze(mean(A1,1));
B2 = logical(A2);
C2 = squeeze(sum(B2,1))+1;
temp2 = squeeze(mean(A2,1))./C2;
temp2a = squeeze(mean(A2,1));
%%
figure;
imshow(uint8(temp1a*13*5))
%
figure;
imshow(uint8(temp2a*33*10))