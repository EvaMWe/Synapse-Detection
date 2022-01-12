%this funtions is a marker-controlled watershed segmentation within
%pre-segmented sections from a whole image containing fluorescence labelled
%synapses
%
% input:
% centroids: indices of centroids wihtin section; 
% section: one pre-segmented section (blobb) containing several connected synapses
% output:
% ListOfAll: 1st col - indices; 2nd cl - related region (integer number)
% S is a threshold to discard dark pixels included by dilatation step
function [ListOfAll,nCinit] = getSingles(centroids, section, S)

section(section < S) = 0;
ListOfAll = find(section ~= 0);
ListOfAll = [ListOfAll zeros(length(ListOfAll),1) zeros(length(ListOfAll),1)];
nC = length(centroids);
nCinit = nC;
%markers = 1:nC;
W = nC +1; %Abgrenzung

%% initialize the dynamic lists
dynList = zeros(length(ListOfAll),2,nC); % pixelID,pixel_value,numberOfCentroids
for c = 1:nC
    %create the initial dynamic list, fill it with centroids neighbors,
    %centroids are labelled with their number
    ListOfAll((ListOfAll(:,1)==centroids(c)),2) = c;
    neighbors = getNeighbors(section,centroids(c),1);
    dynList(1:length(neighbors),1,c) = neighbors;  %pixelIDList
    dynList(1:length(neighbors),2,c) = section(neighbors); %pixelValue
    dynList(:,:,c) = sortrows(dynList(:,:,c),2,'descend');
end

%% start the algorithm
counter = 1:nC;
while ~all(ListOfAll(:,2))
    left_old = length(ListOfAll) - nnz(ListOfAll(:,2));
    % sprintf('%i left',left)
    for c = 1:nC
        clearvars n_old n n_check n_temp
        if c > length(counter)
            continue
        end
        c_ = counter(c);
        point = dynList(1,1,c_);%point with highest intensity in list
        %select list from current dynList for current centroid and sort
        n_old = sort(dynList(1:nnz(dynList(:,1,c_)),1,c_));
        
        n_old(n_old == point) = []; %current point is deleted from List
        if isempty(n_old) %no connected values, one part was separated by thresholding
            ListOfAll(ListOfAll(:,1) == point,2) = W;
            %ListOfAll(ListOfAll(:,1) == point,3) = c_;
            counter(c) = [];
            nC = nC -1;
            if nC == 0
                ListOfAll(ListOfAll(:,2) == 0,2) = W;
            end
            continue
        end
        n = getNeighbors(section,point,1);
        
        
        %check if there are borderpixel
        n_check = n;
        borders = ~ismember(n_check,ListOfAll(:,1));
        if sum(borders) ~= 0
            n(borders) = [];
            ListOfAll((ListOfAll(:,1) == point),2) = W;
            %ListOfAll((ListOfAll(:,1) == point),3) = c_;
        else
            % ckeck the neighbors
            marks = ListOfAll(ismember(ListOfAll(:,1),n_check),2);
            %otherMarks = markers(markers ~= c);
            label = sum(marks ~= c_ & marks ~= 0 & marks~=W);
            if label == 0
                ListOfAll((ListOfAll(:,1) == point),2) = c_;
            else
                ListOfAll((ListOfAll(:,1) == point),2) = W;
                ListOfAll((ListOfAll(:,1) == point),3) = c_;
            end
        end
        
        %refresh dynamic List
        
        % clean list of new points
        replicates = ~ismember(n, n_old);
        n = n(replicates);    %just those that are not in n_old        
        already = ListOfAll((ListOfAll(:,2) ~= 0),1);
        n(ismember(n, already)) = []; %Exclude pixel already marked       
        
        n_temp = unique([n_old; n]); %all neighbors, unique, not already marked --> ready to search next max = point
        dynList(1:length(n_temp),1,c_) = n_temp;  %pixelIDList
        dynList(length(n_temp)+1:end,1,c_) = 0;  %pixelIDList
        dynList(1:length(n_temp),2,c_) = section(n_temp); %pixelValue
        dynList(length(n_temp)+1:end,2,c_) = 0;
        dynList(:,:,c_) = sortrows(dynList(:,:,c_),2,'descend');        
    end
%check if algorithm is still running and number of neighbors decrease
%     left_new = length(ListOfAll) - nnz(ListOfAll(:,2));
%     if left_new == left_old
%         ListOfAll(ListOfAll(:,2)==0,2) = W;
%     end
end
%     replace = ListOfAll((ListOfAll(:,2) == W),3);
%     replace(replace == 0) = [];
%     ListOfAll((ListOfAll(:,2) == W)&(ListOfAll(:,3) ~= 0)&(ListOfAll(:,3) ~= 0),2) = replace;
end


