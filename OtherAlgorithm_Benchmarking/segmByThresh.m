function [region, BW] = segmByThresh(img)

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
region = regionprops(BW, img,'Area','MeanIntensity','PixelIdxList','WeightedCentroid');

% Take only regions which lies in between the specified values:
region = region([region.Area] > minarea & [region.Area] < maxarea);

BW = zeros(size(BW));
for k=1:length(region)
    BW(region(k).PixelIdxList) = 1;
end

% figure
% imagesc(BW)
% title('Binary image after clearing')

end