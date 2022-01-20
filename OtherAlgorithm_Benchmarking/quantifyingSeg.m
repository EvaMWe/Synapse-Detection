% this is a function to quantify the segmentation 
% TP = True Positive: detected puncta has a match on ground truth data
% FP = FAlse Positive: no match on ground truth
% FN = no detection on  (= total number of ground truth puncta - TP)
%Precision: TP/(TP+FP)
%recall: TP/(TP+FN)

%% calculate TP: check for each segmented ROI
function [Precision, recall, Fscore, overlapOfSeg,overlapRel] = quantifyingSeg(roiGT, roiSeg, image)
sizImg = size(image);
nbRoiSeg = length(roiSeg);
nbRoiGT = length(roiGT);
nearest = 3 ; %number of nearest ROIs ('the 3 nearest neighbors')



resultMatrix = cell(nbRoiSeg,7);
for roi = 1:nbRoiSeg
        %get coordinates for corresponding ROI
        roiSeg_coor = [roiSeg(roi).WeightedCentroid];
        t_y = roiSeg_coor(1,2);
        t_x = roiSeg_coor(1,1);
        storeDist = zeros(nbRoiGT,1);
        
        %calculate distance to every GT ROI
        for m = 1:nbRoiGT
            roiGT_coor = [roiGT(m).WeightedCentroid];
            m_x = roiGT_coor(1,1);
            m_y = roiGT_coor(1,2);
            storeDist(m,1) = distanceCent(t_x,m_x,t_y,m_y);
        end
        storeDist_sort = sort(storeDist);
        nearest3 = storeDist_sort(1:nearest);
        %get 3 nearest neighbors
        idxNearest = arrayfun(@(a) getIdxNearest(storeDist,a), nearest3); %number of ROIs, NOT the indexed pixel
        
        %to determine if TP 
        templateRoiSeg = false(sizImg);
        templateRoiSeg(roiSeg(roi).PixelIdxList) = true;
        
        Overlap = zeros(nearest,2);
        for k = 1:nearest
            templateRoiGt = false(sizImg);
            idx = idxNearest(k);
            templateRoiGt(roiGT(idx).PixelIdxList) = true;
            temp = templateRoiSeg & templateRoiGt;
            Overlap(k,1) = sum(sum(temp));
            Overlap(k,2) = idx;
        end
        if sum(Overlap(:,1)) == 0  %then no overlapping           
            resultMatrix{roi,1} = roi;
             resultMatrix(roi,2:end) = num2cell(zeros(1,6));
        else
        [~,I] = sort(Overlap(:,1),'descend');
        Overlap_ = Overlap;
        Overlap =Overlap_(I,:);
        roiAreaGt = length(roiGT(Overlap(1,2)).PixelIdxList);
        roiAreaSeg = length(roiSeg(roi).PixelIdxList);
        areaTogether = roiAreaGt + roiAreaSeg;
        overlap_rel = Overlap(1,1)/(areaTogether - Overlap(1,1));
        % //** include false negative to take not detected ROIs into account
        
        %RoiIdx_over 
        resultMatrix{roi,1} = roi;          %roi_Index_segmented
        resultMatrix{roi,2} = Overlap(1,2); %roi_index_associated reference
        resultMatrix{roi,3} = Overlap(1,1); %overlapping pixel/area
        resultMatrix{roi,4} = Overlap(2:3,2);%index of further neighbored ref. ROIs
        resultMatrix{roi,5} = Overlap(1,1) / roiSeg(roi).Area; %percent overlap
        if resultMatrix{roi,5} > 0.1
            resultMatrix{roi,6} = 1;   %= TP
        else
            resultMatrix{roi,6} = 0; %=FP
        end
        resultMatrix{roi,7} = overlap_rel;
        end
end
resultMatrix = checkFordouble(resultMatrix);
% if idx_doubles ~= 0
%     resultMatrix(idx_doubles,6)= num2cell(2); %ROI is double: currently it will just deleted
% %     overlapping=[cell2mat(resultMatrix(idx_doubles,3)) idx_doubles'];
% %     overlapping = sortrows(overlapping,'descend');
% end

%% calculate parameters
resultNb = cell2mat(resultMatrix(:,6)) ;
%overlap = cell2mat(resultMatrix(:,3));
TP = sum(resultNb == 1);
FP = sum(resultNb == 0);
totalSeg = length(resultNb);
totalGT = length(roiGT);
FN = totalGT - TP;

%//** not detected ROIs are included in relative overlap (overlap_rel)

overlap_extended = [cell2mat(resultMatrix(:,7));zeros(FN,1)];

        


% Precision: TP/(TP+FP)%recall: TP/(TP+FN)
Precision = TP/(TP+FP);
recall = TP/(TP+FN);
Fscore = (2*Precision*recall)/(Precision + recall);
overlapOfSeg = mean(cell2mat(resultMatrix(:,5)));
overlapRel = mean(overlap_extended);


end


    function dist = distanceCent(x2,x1,y2,y1)
        dist = sqrt(((x2-x1)^2+(y2-y1)^2));
    end

    function idx = getIdxNearest(liste, value)
        idx= find(liste == value);
    end
    
    function values_cleared = checkFordouble(values)
    vector = cell2mat(values(:,2)); %associated regions from GT
    area = cell2mat(values(:,3));   %overlapping area
    % deal with zeros:
    Idx_Zeros = vector == 0; %zero Elements
    Idx_no0 = find(vector ~= 0);
    NbZeros = sum(Idx_Zeros);   
    
  
    vector_no0 = vector(~(Idx_Zeros)); %without zeros
    vector_uni = unique(vector_no0);  % unique vlaues, this is sorted
    
    if length(vector_uni) ~= length(vector_no0) %without 0;
        indivDouble = arrayfun(@(a) sum(vector_no0 == a),vector_uni);
        doubleValues = vector_uni(indivDouble > 1);
        nbDoubles = length(doubleValues);
        
        
        %% create new data cell
        values_cleared = cell(length(vector_uni)+NbZeros,size(values,2));
        
        %fill new cell with unique values:
        uniVec = find(vector ~= 0 & ~ismember(vector,doubleValues));
        values_cleared(1:length(uniVec),:) = values(uniVec,:);
        
        %fill new cell with doubles:
        for d = 1:nbDoubles
            x = doubleValues(d);
            idx_x = find(ismember(vector,x));
            values_x = area(idx_x);
            [~, idx_remainer] = max(values_x);
            idx_remainer = idx_x(idx_remainer);
            values_cleared(length(uniVec)+d,:) = values(idx_remainer,:);
            
        end
        
        %fill new cell with Zeros
        values_cleared(length(uniVec)+nbDoubles+1 : end,:) = num2cell(0);
    else
        values_cleared = values;
    end
    end
