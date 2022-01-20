%this function compares the segmenetation quality measured by jaccard and
%bf score
%BW_stack: layer 1 .--> thresholding
%           layer2 --> water shed
%           layer 3 --> local maxima

function [simJac, simBF, Precision, recall, Fscore, overlapOfSeg,overlapRel] = evalSegmentation_woCrop(croppedI,L_ref,BW_ref)
border = 5; %check if border is is equal to border in cellDetection_combined
nbMeth = 4;
simJac = zeros(nbMeth,1);
simBF = zeros(nbMeth,1);
simDice = zeros(nbMeth,1);
Precision = zeros(nbMeth,1);
recall = zeros(nbMeth,1);
Fscore = zeros(nbMeth,1);
overlapOfSeg =zeros(nbMeth,1); 
overlapRel = zeros(nbMeth,1);

simJac_label = zeros(nbMeth,1);

BW_stack = zeros(size(croppedI,1),size(croppedI,2),nbMeth);
L_stack = zeros(size(croppedI,1),size(croppedI,2),nbMeth);
region_stack = cell(nbMeth,1);


%% ground truth
%create BW = groundTruth
%[BW,~,L]  = segmentImage(croppedI);
reg_ref = regionprops(L_ref,croppedI, 'PixelIdxList', 'PixelList', 'Perimeter',...
 'WeightedCentroid','MeanIntensity', 'MaxIntensity');

%% Method 1: Segmentation by Threshold (automated)
[regionsThresh,BW_segmByThresh] = segmByThresh(croppedI);
% set broder to 0
BW_segmByThresh([1:border end-border:end],:) = 0;
BW_segmByThresh(:,[1:border end-border:end]) = 0;

BW_stack(:,:,1) = BW_segmByThresh;
L_stack(:,:,1) = bwlabel(BW_stack(:,:,1));
region_stack{1,1} = regionsThresh;

%% Method 2: Segmentation by Threshold with subsequent water shed separation
[regionsThreshWs,BW_segmByThresh_ws] = segmByThresh_ws(croppedI);
% set broder to 0
BW_segmByThresh_ws([1:border end-border:end],:) = 0;
BW_segmByThresh_ws(:,[1:border end-border:end]) = 0;
BW_stack(:,:,2) = BW_segmByThresh_ws;
region_stack{2,1} = regionsThreshWs;

blank = zeros(size(croppedI));
for roi = 1:length(regionsThreshWs)
    pxlist = regionsThreshWs(roi).PixelIdxList;
    blank(pxlist) = roi;
end
L_stack(:,:,2) = blank;

%% Method 3: Segmentation by local maxima according to Sbalzerinis feature detection 
[regionsLocMax,~,BW_featureDetectionSb] = featureDetectionSb(croppedI, 3, 5, 3, 0);
% set broder to 0
BW_featureDetectionSb([1:border end-border:end],:) = 0;
BW_featureDetectionSb(:,[1:border end-border:end]) = 0;
BW_stack(:,:,3) = BW_featureDetectionSb;
region_stack{3,1} = regionsLocMax;
L_stack(:,:,3) = bwlabel(BW_stack(:,:,3));

%% Method 4: Segmentation by edge filter and advanced post processing
[regionEdgeWs,BW_stack(:,:,4),blank] = cellDetection_iterative_v1(croppedI, 6, 200, 5);
region_stack{4,1} = regionEdgeWs;
%blank = zeros(size(croppedI));
% for roi = 1:length(regionEdgeWs)
%     pxlist = regionEdgeWs(roi).PixelIdxList;
%     blank(pxlist) = roi;
% end
L_stack(:,:,4) = blank;

%% Quantification
for m = 1:nbMeth
    simJac(m,1) = jaccard(BW_ref,logical(BW_stack(:,:,m)));
    simJac_label(m,1) = mean(jaccard(L_ref,L_stack(:,:,m))); % label not properly associated between images (?)
    simBF(m,1) = bfscore(BW_ref,logical(BW_stack(:,:,m))); 

    simDice(m,1) = dice(BW_ref,logical(BW_stack(:,:,m)));
    
    [Precision(m,1),recall(m,1),Fscore(m,1),overlapOfSeg(m,1), overlapRel(m,1)]...
    = quantifyingSeg(reg_ref, region_stack{m,1}, croppedI);
end
end

% % % 
% Lrgb1 = label2rgb(L_stack(:,:,2),'jet','w','shuffle');
% imshow(Lrgb1)
% title('Colored Watershed Label Matrix')


