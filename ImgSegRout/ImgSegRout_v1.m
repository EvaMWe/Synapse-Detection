
% 12.01.2022: release beta v1
% imgSegRout: a function to read out fluorescence traces from image stacks in time lapse 
% microscopy automatically
% Segmentation is performed using SynEdgeWs
% 

function ImgSegRout_v1(preprocessing, iteration, sizeSettings, FrameRange, Func)

alignCommand = Func.Align;
bleachingCorr = Func.bleaching;

[minArea,maxArea] = convert2pixel(sizeSettings);
 
discard = preprocessing.Discard;
denoise = preprocessing.Denoise;
showImage = preprocessing.ShowImage;

it = iteration.it;
%% (1) Data Input and save path


% a) Select experiments
%---------------------------------
    [dataFiles, dataPath] = uigetfile ('.tif','select experiment files', 'Multiselect', 'on'); 
    dataComplete = fullfile(dataPath,dataFiles);
   

% d) name Excel-Sheet
[reportFile,reportPath] = uiputfile('.xlsx','Name Report File');


if ~iscell(dataComplete)
    dataComplete = cellstr(dataComplete);
end

N = length(dataComplete); 


%%
%preallocating data sets
dataSet = cell(N,6);
for exp = 1:N    
    if ~iscell(dataFiles)
        dataFiles = cellstr(dataFiles);
    end
    name = dataFiles{1,exp}; 
    name = replace(name,'.tif','');
    cyStack = LoadMultipage(dataComplete{exp}, 0);  %load image stack
    if alignCommand  == 1
        cyStack = aligneSequence(cyStack);
         sprintf('alignment of experiment %i completed',exp)
    end
   
    %% (2) Find ROIs;
    % a) prepare image for searching ROI,take the last frames after
    % ammonium chloride
    len = size(cyStack,3);
    %range = len-4:len;
    range = FrameRange.Frame1 : FrameRange.Frame2;
    averageImg = mean(cyStack(:,:,range),3,'native');
    % restore averageImg for more precise detection
    averageImgPro = PreProcessing(averageImg,discard,denoise,showImage);
    [regionProp, ~] = cellDetection_iterative_v1_ImgSegRout(averageImgPro, minArea, maxArea, it);
    
    %% (3)Readout and calculation of background traces
    regNb = length(regionProp);
    %restoration of traces //Normalization and Backgroundsubtraction
    [data, backgroundTrace, regionNb] = Readout(regionProp, regNb, cyStack);

    %% (4) Subtraction of background
    %regionNb = size(data,1);
    data_BGsubt = zeros(size(data));
    for region = 1:regionNb
        traceTemp = data(region,:) - backgroundTrace;
        data_BGsubt(region,:) = traceTemp;
    end      
    mean_BG = mean(data_BGsubt,1);
       
    %% (5) smooth traces
    tracesSmth = zeros(size(data));
    regNb = regionNb;
    for i = 1:regNb
        tracesSmth(i,:) = smooth(data_BGsubt(i,:));
    end
    meanCurve = mean(tracesSmth,1);
    
    %% (6) correction for bleaching
    
    % get coefficient
    if bleachingCorr == 1
        [lambda, ~] = createBleachingCurve_woPoly(meanCurve, 13 , 4, 0, 1);  %original, range, cuttingWin, nStim, show images
        
        % perform correction
        dataCorr = zeros(size(data));
        dataCorrSmth = zeros(size(data));
        for k = 1:regNb
            dataCorr(k,:) = itDeconv (data_BGsubt(k,:), lambda);
        end
        for k = 1:regNb
            dataCorrSmth(k,:) = itDeconv (tracesSmth(k,:),lambda);
        end
        
    else
        dataCorr = data_BGsubt;
        dataCorrSmth = tracesSmth;
    end
    
    %meanCorr_nc = mean(dataCorrSmth,1);
    %meanCorr = mean(dataCorrSmth,1); 
      
   sprintf('preprocessing of measurement%i succeeded',exp)
    
        
    
 %----------DATA STORAGE---------------------------------------            
        %% data sheets: rawdata, backgroundsubtracted Data, + smoothed, + bleching corrected, - smoothed)
        datasheet = cat(3,data,data_BGsubt, tracesSmth, dataCorrSmth, dataCorr);  
            
           
    
        
        %% (12) save results as excel file // one per experiment
        filename = replace(dataFiles{exp},'.tif',sprintf('_%s',reportFile));
        export2Excel_afterSeg_Imaging(datasheet,reportPath,filename)
        varnam = replace(dataFiles{exp},'.tif','');
        dataSet(exp,1) = cellstr(varnam);
        dataSet{exp,2} = data;
        dataSet{exp,3} = data_BGsubt;
        dataSet{exp,4} = tracesSmth;
        dataSet{exp,5} = dataCorrSmth; %if bleaching == 0 dataCorrSmth = tracesSmth
        dataSet{exp,6} = dataCorr; %if bleaching == 0 dataCorr = tracesSmth
end

   matname = replace(reportFile,'.xlsx','matfile_');
   cd (reportPath)
   save(matname,'dataSet')
   
end

