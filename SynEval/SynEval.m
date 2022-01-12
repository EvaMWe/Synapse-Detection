%SynEval 
%Version 1.0


function [resultCell] = SynEval_v1(varargin)

%% optional specification of pixel area

nbCostainings = 2;
%% load all image; 
[SynInfo, SynPath] = uigetfile('*.tif','select image for syn-staining (1)',...
    'MultiSelect','on');
syn_Complete = fullfile(SynPath,SynInfo);

[CoStainInfo1, CostainPath1] = uigetfile('*.tif','select images for staining 2',...
    'MultiSelect','on');
    Co1_Complete = fullfile(CostainPath1,CoStainInfo1);

if nbCostainings == 2
     [CoStainInfo2, CostainPath2] = uigetfile('*.tif','select images for staining 3',...
        'MultiSelect','on');
    Co2_Complete = fullfile(CostainPath2,CoStainInfo2);
end

if ~iscell(Co2_Complete)
    Co2_Complete = cellstr(Co2_Complete);
end


if ~iscell(Co1_Complete)
    Co1_Complete = cellstr(Co1_Complete);
end

if ~iscell(Co1_Complete)
    Co1_Complete = cellstr(Co1_Complete);
end

if ~iscell(syn_Complete)
   syn_Complete = cellstr(syn_Complete);
end

imageComplete = cat(1,Co1_Complete,Co2_Complete);


round = size(imageComplete,2); %number of images per staining

resultCell = cell(round+1,8);
featureNames = {'total#Syn', 'fraction Syt1+/Syn', 'FI Syt1+(Syn+)','total2nd', '2nd/Syn', 'FI 2nd(Syn+)',...
    'Syt1/Syn&2nd', 'FI(syt(Syn+2nd)'};
resultCell(1,:) = featureNames;

for r = 1:round
    %% preprocessing syn / other maincostaining
    if ~iscell(syn_Complete)
        syn_Complete = cellstr(syn_Complete);
    end
     imgSyn = double(imread(syn_Complete{1,r}));
     syn_pro = PreProcessing(imgSyn,0,0,0);
     
     %% preprocessing co-stainings
    imgStackCo =zeros(size(imgSyn,1), size(imgSyn,2), nbCostainings );  
    for int = 1:nbCostainings       
        img = double(imread(imageComplete{int,r})); 
        imgStackCo(:,:,int) = PreProcessing(img,0,0,0);
    end

    

%% create template based on syn    
    
 %sprintf('%i',r)
[regionList, binData] = cellDetection_iterative(syn_pro);
    
    %% get thresholds
    threshold_co = zeros(nbCostainings,1);
    for img = 1:nbCostainings
        threshold_co(img,1) = getthreshold(imgStackCo(:,:,img));
    end
    
    %% check regions 
    resultMatrix = zeros(length(regionList),4); 
    % first row: number of region
    % 2° row: mean intensity on syt1 or other...
    % 3t row:mean intensity on vGAT/vGLUT or other... 
    % 4th row:positive for syt.1 and second marker
    % regions below threshold are indicated as zeros
    for img = 1:nbCostainings
        thresh = threshold_co(img,1);
        I = imgStackCo_ori(:,:,img);
        for reg = 1:length(regionList)
           value = mean(mean(I(regionList(reg).PixelIdxList)));
           if value < thresh
               value=0;
           end           
           resultMatrix(reg,img+1) = value;
        end
    end
    doublePositives = (resultMatrix(:,2) ~= 0) & (resultMatrix(:,3) ~= 0);
    resultMatrix(:,4)= resultMatrix(:,2);
    resultMatrix(:,4)=resultMatrix(:,4).*doublePositives;
    syn_nb = length(resultMatrix);
    second_nb = sum(resultMatrix(:,3) ~= 0);  
       
    resultCell{r+1,1} = syn_nb;
    resultCell{r+1,2} = (sum(resultMatrix(:,2) ~= 0))/syn_nb; 
    resultCell{r+1,3} = mean(resultMatrix((resultMatrix(:,2)~= 0),2));
    resultCell{r+1,4} = second_nb; 
    resultCell{r+1,5} = second_nb/syn_nb;
    resultCell{r+1,6} = mean(resultMatrix((resultMatrix(:,3)~= 0),2));
    resultCell{r+1,7} = (sum(resultMatrix(:,4) ~= 0))/second_nb;
    resultCell{r+1,8} = mean(resultMatrix((resultMatrix(:,4)~= 0),2));
end  

end



    
    


