%This is the simple threshold based segmentation improved with a following
%water shed allgorithm...

function [regionStats, BW_ws_cleared] = ThrWs(img)

minarea = 6;           
maxarea = 200;         
pth = 20;              % --> Percentile für die Bestimmung des thresholds



% Determining upper pth-th percentile of intensity values
pth = pth / 100;
[counts, bins] = imhist(img);
k = 0;
while sum(counts(end-k:end)) < pth*sum(counts)
    k = k + 1;
end
level = bins(end-k+1);

% Schwarz-weiß Bild mit 'level' als Schwelle bestimmen:
BW = imbinarize(img, level);

BW = imfill(BW, 'holes');

% Remove regions, that are in contact with border:
BW = imclearborder(BW);

% figure
% imagesc(BW)
% title('Binary image after thresholding')

% apply conditions to refine selection
region = regionprops(BW, img,'Area','MeanIntensity','PixelIdxList','PixelList');

% Take only regions which lies in between the specified values:
region = region([region.Area] > minarea & [region.Area] < maxarea);

BW = zeros(size(BW));
for k=1:length(region)
    BW(region(k).PixelIdxList) = 1;
    
end

%% since here threshold based together with watershed
RegPropNew = region;
nRoi = length(region);
N = length(nRoi)*10;
regionList = repmat(struct('PixelIdxList',[]),1,N);
counter = 1;
unchanged = 0;
for roi = 1:nRoi
    %
    %sprintf('%i _start',roi)
    if RegPropNew(roi).Area >= 40
        region = RegPropNew(roi).PixelList;
        idx = RegPropNew(roi).PixelIdxList;
        [centroids, section, xmin, ymin] = getCentroids(region,idx,img);
        [idxList, numbC] = getSingles(centroids, section,level);
        if length(centroids) > 100
            continue
        end
        idxListtrans = idxList;
        idxList = coordTrans(idxListtrans, size(section),size(BW), xmin, ymin);
        [counter, regionList] = precising(regionList,idxList,numbC, size(BW), counter);
    else
        % condition above threshold:
        regionList(counter).PixelIdxList = RegPropNew(roi).PixelIdxList;
        regionList(counter).Area = length(RegPropNew(roi).PixelList);
        regionList(counter).PixelList = RegPropNew(roi).PixelList;
        counter = counter + 1;
        unchanged = unchanged + 1;
        regionList(counter-1).unchanged = 1;
    end
    
end
regionList = regionList(1:counter-1);

BW_watershed = zeros(size(BW));

%nRoi = length(RegPropNew);
nRoi = length(regionList);
for roi = 1:nRoi
    BW_watershed (regionList(roi).PixelIdxList) = 1;
end
BW_watershed = imbinarize(BW_watershed);
BW_watershed = imfill(BW_watershed,4,'holes');

%convert to greyscale;

regionStats = regionprops(BW_watershed,img, 'PixelIdxList', 'PixelList', 'Area','Perimeter','WeightedCentroid'...
    ,'MeanIntensity','MaxIntensity');
%conditions
regionStats([regionStats.Area] < 6) = [];
%  regionStats([regionStats.MeanIntensity] < 1.0 * thresh) = [];
%  regionStats([regionStats.MaxIntensity]< 1.5 * thresh) = [];

BW_ws_cleared = zeros(size(BW));
nRoi = length(regionStats);
for roi = 1:nRoi
    BW_ws_cleared (regionStats(roi).PixelIdxList) = 1;
end

% figure
% imagesc(BW)
% title('Binary image after clearing')

end