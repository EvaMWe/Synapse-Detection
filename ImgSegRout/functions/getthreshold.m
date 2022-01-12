
%original is not the raw image, it's background subtracted and filtered for
%noise

function [h] = getthreshold(image,original,SN)

mini = min(image(:));
maxi = max(image(:));
image_norm =  (image-mini)/(maxi-mini);


[~,threshold] = edge(image_norm,'sobel');
alpha = -0.2*SN+1.4;
threshold = threshold * alpha;
BW = edge(image_norm,'sobel',threshold);

%remove BW borders
BW([1:11 end-11:end],:) = 0;
BW(:,[1:11 end-11:end]) = 0;
% 
% figure ('Name','after Smoothing and Normalization');
% imshow(image_norm);
%
%visualize detection
% vis_cells(filtered,BW,1,sprintf('after edge filtering for,%s', name),col);

%dilate
SE = strel('disk',2,0);
BWdil = imdilate(BW,SE);

%fill holes
BWfill = imfill(BWdil,4,'holes');

%vis_cells(image_norm,BWfill);

%convert to greyscale;
regionStats = regionprops(BWfill,original,'PixelIdxList','Area');
regionStats([regionStats.Area] < 6) = [];
regionStats([regionStats.Area] > 200 ) = [];

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

