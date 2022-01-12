% this function further calculates parameters for syt1 uptake experiments
% from MergeMarker within Syt1Eval

function [totalSyt_pos, fraction_pos,totalSyt_pos2, fraction_pos2,QualityFraction,IntSyt_pos,IntSyt_pos2]  = postAnalysis(results,results_control, totalNumb)

totalSyt_ = size(results,1)-1; %not validated by another marker
totalM1 = totalNumb(1);
totalM2 = sum(cell2mat(results_control(2:end,2)));
M1 = cell2mat(results(2:end,:,1)); %array for marker 1
M2 = cell2mat(results(2:end,:,2)); %array for marker 2

logicalM1 = logical(M1(:,2));
IntensityM1 = M1(:,end);
logicalM2 = logical(M2(:,3));
IntensityM2 = M2(:,end);


totalSyt_pos = sum(logicalM1); %Ppositiv for Syn;
totalSyt_pos2 = sum(logicalM1 & logicalM2); %positive for both marker

IntSyt_pos = mean(IntensityM1(logicalM1));
IntSyt_pos2 = mean(IntensityM2(logicalM1 & logicalM2));%positive for both marker


fraction_pos = totalSyt_pos/totalM1;
fraction_pos2 = totalSyt_pos2/totalM2;

QualityFraction = totalSyt_pos / totalSyt_;


end



