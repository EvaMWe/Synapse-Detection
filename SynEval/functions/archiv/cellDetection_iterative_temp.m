%% Segmentation of synaptic puncta by iterative thresholding
%last changes/update: 06.01.02022: merge regionList after each iteration
function [regionStats,BW_watershed, L_watershed] = cellDetection_iterative_temp(image,minArea,maxArea,it)
close all
%% get label of figures
%default
% name = 'unspecified';
% col = 'red';
border = 5; %should be broader than kernel width in PreProcessing_imaging

% if nargin >= 5
%     name = varargin{2};
% end
% 
% if nargin == 6
%     col = varargin{3};
% end


%% just if it's not yet normalised; won't change anything if it already is
mini = min(image(:));
maxi = max(image(:));
filtered =  (image-mini)/(maxi-mini);

%% calculate threshold
% ImgVec = filtered(:);
% ImgVec = double(ImgVec);
% ImgVec(ImgVec==0) = [];
% [counts,bins] = hist(ImgVec,500);
% %hist(ImgVec,50) falls mans doch mal sehen will
% k = 0;
% pth = 0.7;
% while sum(counts(end-k:end))< pth*sum(counts) %threshold finden
%     k = k+1;
% end
% thresh=bins(end-k+1);

%%

s = it;
[~,~,Gv,Gh] = edge(filtered,'sobel');
%RegionStorage = cell(s,2);
Op =sqrt((Gv.*Gv)+(Gh.*Gh));
hInit =max(Op(:));

pth_0 = 0.9;
thresh = calcthresh(filtered,pth_0);
for step = 1:s
    x = step;
    pth = 0.01*x;    
    [BW, h] = S(Op,pth, hInit);
    %remove BW borders
    BW([1:border end-border:end],:) = 0;
    BW(:,[1:border end-border:end]) = 0;
    
   
    % figure ('Name','after Smoothing and Normalization');
    % imshow(filtered);
    % %
    %visualize detection
    % vis_cells(filtered,BW,col,1,sprintf('after edge filtering for,%s', name));
    
    %dilate
    SE = strel('disk',1,0);
    BWdil = imdilate(BW,SE);
    BWfill = imfill(BWdil,4,'holes');
    BWer = imerode(BWfill,SE);
    %fill holes
    BWfill = BWer;
    BWfill = logical(BWfill);
    %     Nozero = filtered(BWfill);
    %     Nozero(Nozero == 0) = [];
    %     h = Minimum_median(Nozero, 5, 'Type','percent', 'Dimension',1 );
    hInit = h;
    
    
    %vis_cells(filtered,BWfill,col,1,sprintf('after dilatation and flood filling,%s', name));
    
    %get regiones
    
    RegProp = regionprops(BWfill, 'PixelIdxList', 'PixelList', 'Area','Perimeter');
    
    % delete small regions;
    RegPropNew = RegProp;
    RegPropNew([RegProp.Area] <= floor(minArea*0.5)) = [];
    
    BWcleared = zeros(size(BW));
    nRoi = length(RegPropNew);
    for roi = 1:nRoi
        BWcleared(RegPropNew(roi).PixelIdxList) = 1;
    end
    
    BWcleared = imbinarize(BWcleared);
    %vis_cells(filtered,BWcleared,col,1,sprintf('after clearing,%s', name));
    
    %particle refinement
    N = length(RegProp)*10*s;
    regionList = repmat(struct('PixelIdxList',[]),1,N);
    counter = 1;
    unchanged = 0;
    
    for roi = 1:nRoi
        %
        %sprintf('%i',roi)
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
            
            [idxList, numbC] = getSingles(centroids, section,thresh);
            
            idxListtrans = idxList;
            idxList = coordTrans(idxListtrans, size(section),size(BWfill), xmin, ymin);
            [counter, regionList] = precising(regionList,idxList,numbC, size(BWfill), counter,h, filtered, minArea);
            % BWcur = zeros(size(BW));
            %  BWcur(regionList(counter-3).PixelIdxList) = 1;
            % vis_cells(image_norm ,BWcur,col,1,sprintf('after watershed,%s', name));
        else
            %        condition above threshold:
            PixelIdxList = RegPropNew(roi).PixelIdxList;
            regionList(counter).Area = length(RegPropNew(roi).PixelList);
            regionList(counter).PixelList = RegPropNew(roi).PixelList;
            
            % discard "tails'
            BWcur = zeros(size(BW));
            BWcur(PixelIdxList) = 1;
            if length(PixelIdxList) >= minArea*2
                BWdisc = imerode(BWcur,SE);
                BWdisc = imdilate(BWdisc,SE);
            else
                BWdisc = BWcur;
            end
            regiondisc = regionprops(BWdisc,'Area','PixelList','PixelIdxList');
            %sprintf('%i',roi)
            if size(regiondisc,1) >= 6
                %sprintf('%i',roi)
                [regionList, counter] = checkIfMultiple(regionList, regiondisc, counter, BW);
            elseif size(regiondisc,1) == 0
                regionList(counter) = [];
                counter = counter + 1;
            else
                regionList(counter).PixelIdxList = regiondisc.PixelIdxList;
                counter = counter + 1;
            end
            
            %
            unchanged = unchanged + 1;
            regionList(counter-1).unchanged = 1;
            
            
            %vis_cells(filtered ,BWdisc);
        end
        
    end
    
    %initial total list
    
    regionList = regionList(1:counter-1);
    if size(regionList,2) == 0
        continue
    end
    regionList([regionList.Area] < minArea) = [];
    regionList(cellfun(@isempty,{regionList.Area})) = [];
    regionList(cellfun(@isempty,{regionList.PixelIdxList})) = [];
    %%current BW template
    BW_watershed = zeros(size(BW));
    nRoi = length(regionList);
    for roi = 1:nRoi
        BW_watershed (regionList(roi).PixelIdxList) = 1;
    end
    
    
    %% merge lists
    
    if step == 1
        nReg = length(regionList);
        regionList_pre =  zeros(maxArea*2,nReg*s*2);
        for roi = 1:nReg
            regionList_pre(1:length(regionList(roi).PixelIdxList),roi) = regionList(roi).PixelIdxList;
        end
       
        %regionList_pre(cellfun(@isempty,{regionList_pre.PixelIdxList})) = [];
    else
        
        nReg = length(regionList); %number of new regions
       
        for roi = 1:nReg
            %write current list to temporary storage            
            %fill with regionList_pre
            RegpreFilled = regionList_pre;
            RegpreFilled( :, ~any(RegpreFilled,1) ) = []; %remove columns that are not yet filled
            nRegpre = size(RegpreFilled,2);
                       
            %create current label matrix
            L_pre= zeros(size(BWcleared));
            numPixel = zeros(nRegpre,1);
            
            for ro = 1:nRegpre
                RegpreFilledCol = RegpreFilled(:,ro);
                pxlList = RegpreFilledCol(RegpreFilledCol~=0);
                L_pre(pxlList) = ro;
                numPixel(ro,1) = length(pxlList);
            end
            
            reg = regionList(roi).PixelIdxList;
            regionNr = unique(L_pre(reg)); %label in L_pre is associated to current regionList_pre
            regionNr = regionNr(regionNr ~= 0); %regions that are overlapping with current region
            if ~isempty(regionNr)
                %refine regions
                partNb = length(regionNr);
                partPixlList = zeros(200,partNb+1);
                for part = 1:partNb+1
                    if part == 1
                        partPixlList(1:length(reg),part) = reg;
                    else
                        regTemp = RegpreFilled(:,regionNr(part-1));
                        partPixlList(1:length(regTemp),part) = regTemp;
                    end
                end
                
                [newRegions] = checkRegion(partPixlList, filtered,minArea);
                nbRegions = size(newRegions,2);
                
                %write regions into temporary list
                
                %remove regions in list of number of new lists is below
                %number of orgignally participating lists
                if nbRegions < length(regionNr)
                                           
                    d = (length(regionNr)) - nbRegions; 
                    idx = regionNr((end-d):end);
                    regionList_pre(:,idx) = [];                    
                end
                                
                for r = 1:nbRegions
                   
                    if r > length(regionNr)
                        region = newRegions(:,r);                         
                        regionList_pre(1:length(region),size(RegpreFilled,2)+(r-length(regionNr))) = region; 
                        regionList_pre(length(region)+1:end,size(RegpreFilled,2)+(r-length(regionNr))) = 0;
                         
                        
                    else
                        region = newRegions(:,r);                         
                        regionList_pre(1:length(region),regionNr(r)) = region;
                        regionList_pre(length(region)+1:end,regionNr(r)) = 0;
                       
                    end
                    
                end
                
            else
                %no overlap, simply add region
                regionList_pre(1:length(reg),size(RegpreFilled,2)+1) = reg;
               
                
            end
           
            
        end
        
    end
    
end

% reshape regions;
%remove zero columns:
regionList = regionList_pre;
regionList( :, ~any(regionList,1) ) = [];
nreg = size(regionList,2);


BW_watershed = zeros(size(BW));
L_watershed = zeros(size(BW));
%nRoi = length(RegPropNew);

for r = 1:nreg
    region = regionList(:,r);
    pixelNb = nnz(region);
    region = region(1:pixelNb,1);
    region = cleanRg(region,size(filtered));
    L_watershed (region) = r;
    BW_watershed (region) = 1;
end
regionStats = regionprops(L_watershed,image, 'PixelIdxList', 'PixelList', 'Area','Perimeter','WeightedCentroid'...
    ,'MeanIntensity','MaxIntensity');
% clear
regionStats([regionStats.MeanIntensity] == 0) = [];
regionStats([regionStats.Area] < minArea ) = [];
regionStats([regionStats.Area] > maxArea ) = [];
regionStats(cellfun(@isempty,{regionStats.PixelIdxList})) = [];





%colored watershed
figure
Label = label2rgb(L_watershed,'jet','w','shuffle');
imshow(Label)
title('Colored Watershed Label Matrix')

end

%%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
function [Op2,h] = S(Op, pth, hInit)

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

function [regList, nbCentroids] = checkRegion(partMat,image, minArea)
         fuse = unique(partMat(partMat ~= 0));         
         thresh = min(image(fuse));
         [fusex, fusey] = ind2sub(size(image),fuse);
         fuse_xy = [fusey fusex];
         [centroids, section, xmin, ymin] = getCentroids(fuse_xy,fuse,image);
         nbCentroids = length(centroids);
         if nbCentroids == 1
             regList = fuse;
             
         elseif length(fuse) <= minArea * 4
             regList = fuse;            
             
         else  
             
             [idxList, ~] = getSingles_all(centroids, section,thresh);
             idxListtrans = idxList;
             newRegions = coordTrans(idxListtrans, size(section),size(image), xmin, ymin);
             regList = precising_woList(newRegions,length(centroids), size(image), thresh, image);
         end         
end

function regCleaned = cleanRg(PixelIdxList,sizImg)
SE = strel('disk',1,0);
BWcur = zeros(sizImg);
BWcur(PixelIdxList) = 1;
BWdisc = imerode(BWcur,SE);
BWdisc = imdilate(BWdisc,SE);
BWdisc = imfill(BWdisc,4,'holes');
regCleaned = find(BWdisc ==1);
end

