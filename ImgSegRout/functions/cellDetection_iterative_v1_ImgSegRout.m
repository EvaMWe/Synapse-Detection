% 01.12.2022: release beta v1
% last term of function name declares belonging routine

% this function performs the segmentation of fluorescence images to detect
% synaptic puncta automatically

function [regionStats,BW_watershed, L_watershed] = cellDetection_iterative_v1_ImgSegRout(image,minArea,maxArea,it)
close all
border = 5; %if changed, it should be broader than kernel width in PreProcessing


%% (1) Image Normalization
%just if it's not yet normalised; won't change anything if it already is
mini = min(image(:));
maxi = max(image(:));
filtered =  (image-mini)/(maxi-mini);


%% (2) Initialize segmentation: create gradient image
s = it;
[~,~,Gv,Gh] = edge(filtered,'sobel');
RegionStorage = cell(s,2);
BW_zero = false(size(filtered));
BW_zero([1:border end-border:end],:) = 0;
BW_zero(:,[1:border end-border:end]) = 0;
Op =sqrt((Gv.*Gv)+(Gh.*Gh));

c = 1;
pth_0 = 0.95;
thresh = calcthresh(filtered,pth_0);
hInit =max(Op(:));

%% (3) iterative thresholding followed by morphological operations and watershed transformation
for step = 1:s
    x = step-1;
    pth = 0.005*(1.259)^x;    
    BW = S(Op,pth,hInit);    
    %remove BW borders
    BW([1:border end-border:end],:) = 0;
    BW(:,[1:border end-border:end]) = 0;
    
    BW = BW-BW_zero;
    
    % morphological operators
    SE = strel('disk',1,0);
    BWdil = imdilate(BW,SE);
    BWfill = imfill(BWdil,4,'holes');
    BWer = imerode(BWfill,SE);
    %fill holes
    BWfill = BWer;
    BWfill = logical(BWfill);  
    Nozero = filtered(BWfill);
    Nozero(Nozero == 0) = [];
    h = Minimum_median(Nozero, 5, 'Type','percent', 'Dimension',1 );
    %%
   
    %get regiones    
    RegProp = regionprops(BWfill, 'PixelIdxList', 'PixelList', 'Area','Perimeter');
    
    % delete small regions;
    RegPropNew = RegProp;
    RegPropNew([RegProp.Area] <= minArea) = [];
    
    BWcleared = zeros(size(BW));
    nRoi = length(RegPropNew);
    for roi = 1:nRoi
        BWcleared(RegPropNew(roi).PixelIdxList) = 1;
    end
        
    
    % particle refinement
    N = length(RegProp)*10;
    regionList = repmat(struct('PixelIdxList',[]),1,N);
    counter = 1;
    unchanged = 0;
    for roi = 1:nRoi      
       
        if RegPropNew(roi).Area >= minArea*3
            region = RegPropNew(roi).PixelList;
            idx = RegPropNew(roi).PixelIdxList;
            
            [centroids, section, xmin, ymin] = getCentroids(region,idx,filtered);
            if length(centroids) > 200 || isempty(centroids)
                continue
            elseif length(centroids) == 1
                PixelIdxList = RegPropNew(roi).PixelIdxList;
                % discard "tails'
                BWcur = zeros(size(BW));
                BWcur(PixelIdxList) = 1;
                BWdisc = imerode(BWcur,SE);
                BWdisc = imdilate(BWdisc,SE);
                regiondisc = regionprops(BWdisc,'Area','PixelList','PixelIdxList');
                [regionList, counter] = checkIfMultiple(regionList, regiondisc, counter, BWdisc);
                
                unchanged = unchanged + 1;
                regionList(counter-1).unchanged = 1;
                continue
            end
            
            [idxList, numbC] = getSingles(centroids, section, thresh);
            
            idxListtrans = idxList;
            idxList = coordTrans(idxListtrans, size(section),size(BWfill), xmin, ymin);
            [counter, regionList] = precising(regionList,idxList,numbC, size(BWfill), counter,h, filtered, minArea);
        else
            % condition above threshold:
            PixelIdxList = RegPropNew(roi).PixelIdxList;
            regionList(counter).Area = length(RegPropNew(roi).PixelList);
            regionList(counter).PixelList = RegPropNew(roi).PixelList;
            
            % discard "tails'
            BWcur = zeros(size(BW));
            BWcur(PixelIdxList) = 1;
            if length(PixelIdxList) >= minArea * 2
                BWdisc = imerode(BWcur,SE);
                BWdisc = imdilate(BWdisc,SE);
            else
                BWdisc = BWcur;
            end
            regiondisc = regionprops(BWdisc,'Area','PixelList','PixelIdxList');
           
            if size(regiondisc,1) ~= 0
                [regionList, counter] = checkIfMultiple(regionList, regiondisc, counter, BW);
            else
                regionList(counter) = [];
                counter = counter + 1;
            end
            
            
            unchanged = unchanged + 1;
            regionList(counter-1).unchanged = 1;
            
        end
        
    end
    regionList = regionList(1:counter-1);
    if size(regionList,2) == 0
        continue
    end
    regionList(cellfun(@isempty,{regionList.Area})) = [];
    regionList(cellfun(@isempty,{regionList.PixelIdxList})) = [];
    regionList([regionList.Area] < minArea) = [];
    
    % current BW template
    BW_watershed = zeros(size(BW));    
    nRoi = length(regionList);
    for roi = 1:nRoi        
        BW_watershed (regionList(roi).PixelIdxList) = 1;
    end
    
    BW_zero = logical(BW_zero+BW_watershed);
    
    RegionStorage{c,1} = regionList;
    RegionStorage{c,2} = length(regionList);
    c = c+1;
end


%% (4) create whole list with ROI coordinates
numbRoi = sum(cell2mat(RegionStorage(:,2)));

regionList = repmat(struct('PixelIdxList',[]),numbRoi,1);
count=1;
for slice = 1:c-1    
    regL = RegionStorage{slice,1};
    for roi = 1:length(regL)
        regionList(count+roi-1).PixelIdxList = regL(roi).PixelIdxList;        
    end
    count= count+length(regL)-1;
end

regionList(cellfun(@isempty,{regionList.PixelIdxList})) = [];

regionStats = regionList;

%% (5) clear overlaps and create binary mask (BW_watershed) and labeled matrix (L-watershed)
for clearing = 1:3
    BW_watershed = zeros(size(BW));
    L_watershed = zeros(size(BW));
    %nRoi = length(RegPropNew);
    nRoi = length(regionStats);
    for roi = 1:nRoi
        L_watershed (regionStats(roi).PixelIdxList) = roi;
        BW_watershed (regionStats(roi).PixelIdxList) = 1;
    end
    regionStats = regionprops(L_watershed,image, 'PixelIdxList', 'PixelList', 'Area','Perimeter','WeightedCentroid'...
        ,'MeanIntensity','MaxIntensity');
    % clear
    regionStats([regionStats.MeanIntensity] == 0) = [];
    regionStats([regionStats.Area] < minArea ) = [];
    regionStats([regionStats.Area] > maxArea ) = [];
    regionStats(cellfun(@isempty,{regionStats.PixelIdxList})) = [];
end

%% (6) visualize segmentation
%colored watershed

Label = label2rgb(L_watershed,'jet','w','shuffle');
imshow(Label)
title('Colored Watershed Label Matrix')
vis_cells(filtered ,BW_watershed);


end

%%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function Op2 = S(Op, pth, hInit)

Oplin = Op(:);
[counts,bins] = hist(Oplin,500);
k = 0;
while sum(counts(end-k:end))< pth*sum(counts) %threshold finden
    k = k+1;
end
h=bins(end-k+1);
Op2 = Op;
Op2(Op2<h | Op2>hInit) = 0;
Op2 = logical(Op2);
end

function [regList_all, counter] = checkIfMultiple(regList_all,regList,counter, img) %image just for size
template = zeros(size(img));
counterInit = counter;

if isfield(regList,'PixelIdxList') && ~isempty(regList) 
    template(regList.PixelIdxList)=1;
    select = bwconncomp(template,4);
    if select.NumObjects > 1
        numObj = select.NumObjects;
        for nO = 1:numObj
            PixelIdxList = select.PixelIdxList{1,nO};
            %regCleaned = cleanRg(PixelIdxList,sizImg,SE);
            if length(PixelIdxList) < 6              
                continue
            end
            regList_all(counter).PixelIdxList = PixelIdxList;
            regList_all(counter).Area = length(PixelIdxList);
            [r,c] = ind2sub(size(img),PixelIdxList);
            regList_all(counter).PixelList = [r c];
            counter = counter+1;
        end
        if counter == counterInit % if all regions are below min numb pixel
            counter = counter + 1;
        end
    else        
        regList_all(counter).PixelIdxList = regList.PixelIdxList;
        regList_all(counter).Area = length(regList.PixelIdxList);
        [r,c] = ind2sub(size(img),regList.PixelIdxList);
        regList_all(counter).PixelList = [r c];
        counter = counter+1;
    end
end
end

function S = calcthresh(img,pth)
ImgVec = img(:);
ImgVec = double(ImgVec);
ImgVec(ImgVec==0) = [];
[counts,bins] = hist(ImgVec,500);
%hist(ImgVec,50) falls mans doch mal sehen will
k = 0;
while sum(counts(end-k:end))< pth*sum(counts) %threshold finden
    k = k+1;
end
S = bins(end-k+1);
end

