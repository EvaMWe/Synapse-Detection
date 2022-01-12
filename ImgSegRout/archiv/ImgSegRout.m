%% 31.08.2021 // with watershed algorithm
%% in this a function that will run a routin to evaluate Imaging Recordings with SyPhy or mOrange
% Input: varargin{1} = alignCommand: 
%           0 = do not automatically align image sequence; 
%           1 = automatic alignment is started
%           default = 1;
% varagin{2] = bleaching: 0 = turn correction off; default: 1
% load Data
% Segmentation
% 
function ImgSegRout(varargin)

alignCommand = 0;
bleachingCorr = 0;
%% (1) Data Input and save path

if nargin > 0
    alignCommand = varargin{1};
end

% a) Select experiments
%---------------------------------
    [dataFiles, dataPath] = uigetfile ('.tif','select experiment files', 'Multiselect', 'on'); 
    dataComplete = fullfile(dataPath,dataFiles);
   
% b) some defaults:

N = numel(dataComplete); 


% c) select the directory to save the results
%[SaveName, SavePath] = uiputfile('.mat', 'select directory to save your files');

% d) name Excel-Sheet
[reportFile,reportPath] = uiputfile('.xlsx','Name Report File');


% f) enter bleaching segment
% prompt = {'Enter frame number for bleaching correction'};
% nBleach = inputdlg(prompt);
% nBleach_ = str2double(nBleach{1,1});

if ~iscell(dataComplete)
    dataComplete = cellstr(dataComplete);
end


%%
%preallocating data sets
dataSet = cell(N,6);
for exp = 1:N    
    if ~iscell(dataFiles)
        dataFiles = cellstr(dataFiles);
    end
    name = dataFiles{1,exp}; 
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
    range = 140:142;
    averageImg = mean(cyStack(:,:,range),3,'native');
    % restore averageImg for more precise detection
    averageImgPro = PreProcessing(averageImg,1,1,0);
    [regionProp, ~] = cellDetection_iterative(averageImgPro);
    
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

   matname = replace(reportFile,'.xlsx','MatData');
   cd (reportPath)
   save(matname,'dataSet')
   
end

