

% input is background subtracted image, not normalized

function [h] = getthreshold(image,minArea,maxArea)

original = image;
mini = min(image(:));
maxi = max(image(:));
image_norm =  (image-mini)/(maxi-mini);


[~,threshold] = edge(image_norm,'sobel');
alpha = 0.9;
threshold = threshold * alpha;
BW = edge(image_norm,'sobel',threshold);

%remove BW borders
BW([1:11 end-11:end],:) = 0;
BW(:,[1:11 end-11:end]) = 0;
% 
%dilate
SE = strel('disk',2,0);
BWdil = imdilate(BW,SE);

%fill holes
BWfill = imfill(BWdil,4,'holes');

%vis_cells(image_norm,BWfill);

%convert to greyscale;
regionStats = regionprops(BWfill,original,'PixelIdxList','Area');
regionStats([regionStats.Area] < minArea) = [];
regionStats([regionStats.Area] > maxArea ) = [];

Intensity = zeros(length(regionStats),1);
nRoi = length(regionStats);
for roi = 1:nRoi
    Intensity(roi) =mean(mean(original(regionStats(roi).PixelIdxList)));
end
% Intensity = original(BWfill);
% Intensity(Intensity == 0) = [];
% 
% Intensityrev = original(~BWfill);
% %Intensityrev(Intensityrev == 0) = [];
% 
% h1 = Minimum_median(Intensity, 10, 'Type','percent', 'Dimension',1 );
% h2 = Maximum_median(Intensityrev, 10, 'Type','percent', 'Dimension',1 );
% h3 = (mean(original(:)))* SN * 0.8 ;
  h=Minimum_median(Intensity, 1, 'Type','percent', 'Dimension',1 );
end

