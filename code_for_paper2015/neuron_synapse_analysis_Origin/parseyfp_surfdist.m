path = arg1;

mergedpath = [path 'storm_merged/'];
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];

prepath = [path 'selected_presynaptic/sections/'];
prefiles = [dir([prepath '*.png']) dir([prepath '*.tif'])];

yfppath = [path 'allyfp/'];
%yfppath = [path 'nearconn_yfp/'];
yfpfiles = [dir([yfppath '*.png']) dir([yfppath '*.tif'])];
num_images = numel(yfpfiles);
info = imfinfo([yfppath yfpfiles(1,1).name]);

%set voxel size in nm
voxel=[15.8, 15.8, 70];



   pathout =  [path 'parseyfp/' ];

block = floor(str2num(blocksize)/voxel(1));
overlap = floor(str2num(overlapsize)/voxel(1));
disp(block)
disp(overlap)
disp(num_images)
numxblocks = ceil((info(1,1).Width-2*overlap)/block);
numyblocks = ceil((info(1,1).Height-2*overlap)/block);
disp(numxblocks)
disp(numyblocks)
%
counter = 0;
startendpos = zeros(numxblocks*numyblocks,4);
disp(size(startendpos))
for nx=1:numxblocks
    startxpix = nx*block-block + 1;
    endxpix = nx*block + overlap*2;
    if nx==numxblocks
        endxpix = info(1,1).Width;
    end
    for ny=1:numyblocks
        counter = counter+1;
        startypix = ny*block-block + 1;
        endypix = ny*block + overlap*2;
        if ny==numyblocks
            endypix = info(1,1).Height;
        end
        
        startendpos(counter,:)= [startxpix endxpix startypix endypix];   
    end
    disp(counter)
end

xpix = zeros(numel(startendpos(:,1)),startendpos(1,2)-startendpos(1,1)+1);
ypix = zeros(numel(startendpos(:,1)),startendpos(1,4)-startendpos(1,3)+1);
%
for i=1:numel(startendpos(:,1))
ypix(i,(1:(startendpos(i,4)-startendpos(i,3)+1))) = (startendpos(i,3):startendpos(i,4));
xpix(i,(1:(startendpos(i,2)-startendpos(i,1)+1))) = (startendpos(i,1):startendpos(i,2));
end
    gausspix = 3;%str2num(gauss);
    sigmaZ = gausspix* voxel(1)/voxel(3);
    sizZ= gausspix*2*sigmaZ* voxel(1)/voxel(3);
    xZ=-ceil(sizZ/2):ceil(sizZ/2);
    H1 = exp(-(xZ.^2/(2*sigmaZ^2)));
    H1 = H1/sum(H1(:));
    Hz=reshape(H1,[1 1 length(H1)]);

    i = str2num(slicenum);
    curr_startend = startendpos(i,:);
    %load in file info
    %allocate arrays
    disp('allocating arrays')
    BY = zeros(size(xpix,2), size(ypix,2), num_images,'uint8');

    disp('loading y')
    for j=1:num_images
        disp(j)
        yp=ypix(i,find(ypix(i,:)));
        xp=xpix(i,find(xpix(i,:)));    
        AA = imread([yfppath yfpfiles(j,1).name]);
        AA = AA(yp,:);
        AA = AA(:,xp);
        BY(1:size(AA,1),1:size(AA,2),j)= AA(:,:);
    end
    bsub = [];
    if numel(find(BY))>0
    %by2 = zeros(size(BY));
    %for k=1:num_images
    %    disp(k)
    %   by2(:,:,k) = mat2gray(imfilter(BY(:,:,k),fspecial('gaussian',gausspix*2,gausspix),'same','replicate'));
    %end
    %by2 = imfilter(by2,Hz,'same','replicate');% show blurred images

     %create surface of neuron and distance matrix
    B=imfill(BY>(0.2),'holes');
    B1 = bwperim(B,26);
    D = bwdistsc1(B1,[voxel(1),voxel(2),voxel(3)],2000);
    D(B>0)=-D(B>0);

    
    %these save out B B1 and D in list form for each block
    B = B((overlap+1):(end-overlap),(overlap+1):(end-overlap),:); 
    [bsub2, bsub1, bsub3] = ind2sub(size(B),find(B));
    bsub = cat(2,bsub1,bsub2,bsub3);
    bsub(:,1)= bsub(:,1) + overlap + curr_startend(1) - 1;
    bsub(:,2)= bsub(:,2) + overlap + curr_startend(3) - 1;
    
    B1 = B1((overlap+1):(end-overlap),(overlap+1):(end-overlap),:);     
    [b1sub2, b1sub1, b1sub3] = ind2sub(size(B1),find(B1));
    b1sub = cat(2,b1sub1,b1sub2,b1sub3);
    b1sub(:,1)= b1sub(:,1) + overlap + curr_startend(1) - 1;
    b1sub(:,2)= b1sub(:,2) + overlap + curr_startend(3) - 1;
    parsave([pathout sprintf('B1_%d.mat', i)], b1sub);  

    D = D((overlap+1):(end-overlap),(overlap+1):(end-overlap),:); 
    [dsub2, dsub1, dsub3] = ind2sub(size(D),find(D<2000 & D>-2000));
    dsub4 = zeros(size(dsub1));
    for j=1:numel(dsub1)
        dsub4(j,1) = D(dsub2(j), dsub1(j), dsub3(j));
    end
    dsub = cat(2,dsub1,dsub2,dsub3,dsub4);
    dsub(:,1)= dsub(:,1) + overlap + curr_startend(1) - 1;
    dsub(:,2)= dsub(:,2) + overlap + curr_startend(3) - 1;   
    parsave([pathout sprintf('D_%d.mat', i)], dsub);
    end
    	parsave([pathout sprintf('B_%d.mat', i)], bsub);
