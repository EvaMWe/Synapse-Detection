function export2Excel_afterSeg_staining(data,savePath,filename)
cd (savePath)


id = 'MATLAB:xlswrite:AddSheet';
warning('off',id);

writetable(array2table(data),filename,'WriteVariableNames',false,...
    'Sheet','data');


end
