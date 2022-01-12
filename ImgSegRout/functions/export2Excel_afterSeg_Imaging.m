function export2Excel_afterSeg_Imaging(data,savePath,filename)
cd (savePath)


id = 'MATLAB:xlswrite:AddSheet';
warning('off',id);


writetable(array2table(data(:,:,1)),filename,'WriteVariableNames',false,...
    'Sheet','raw_data');
writetable(array2table(data(:,:,2)),filename,'WriteVariableNames',false,...
    'Sheet','after BGsubtr');
writetable(array2table(data(:,:,3)),filename,'WriteVariableNames',false,...
    'Sheet','smoothed');
writetable(array2table(data(:,:,4)),filename,'WriteVariableNames',false,...
    'Sheet','corrected for bleaching');
writetable(array2table(data(:,:,5)),filename,'WriteVariableNames',false,...
    'Sheet','corrected not smoothed');

end
