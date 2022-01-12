% compute the x highest values in an array
% Version from 18.10.2019
%SYNTAX
%[med_highest,maxs] = Maximum_median(array,amount);
%[med_highest,maxs] = Maximum_median(array,amount,...);
%
%DESCRIPTION
%input:
%array: data
%amount: fractionnumber or percent of highest values
%
% USE NAME VALUE PAIRS 
%* 'Dimension' - n 
%  --> dimension of sorting
%     1 = sorts allong col
%     2 = sorts alling row
%* 'Type' - 'frac', 'percent', 'number' 
%   --> type = define the number of highest values
%       'frac': fraction (number = frac * size(array,DIM), 
%       'percent': percentage (number = percent/100*size(array,DIM)
%       'number': number = number
%

function [med_highest,maxs] = Maximum_median(array,amount,varargin)
%% set default
dim = 1;
type = 'frac';

%% read in input arguments
if nargin > 2
    nArg = length(varargin);    
    for arg = 1:2:nArg
        nameArg = varargin{arg};
        switch nameArg
            case 'Dimension'
                dim = varargin{arg+1};
            case 'Type'
                type = varargin{arg+1};
            otherwise
        end
    end
end
    

if dim == 1  %sorts elements in a column
    if strcmp(type,'frac')
        numb_highest = round(amount*size(array,1));
        if numb_highest == 0
            numb_highest = 1;
        end
    elseif strcmp(type,'percent')
        numb_highest = round(amount/100*size(array,1));
    elseif strcmp(type,'number')
        numb_highest = amount;
    else
        error('invalid input argument')
    end
    array_sort = sort(array,1, 'descend');
    maxs = array_sort(1:numb_highest,:);
    med_highest(1,:) = median(maxs,1);
    
    
elseif dim == 2  %sorts elements in a raw
    if strcmp(type,'frac')
        numb_highest = round(amount*size(array,2));
        if numb_highest == 0
            numb_highest = 1;
        end
    elseif strcmp(type,'percent')
        numb_highest = round(amount/100*size(array,2));
    elseif strcmp(type,'number')
        numb_highest = amount;
    else
        error('invalid input argument')
    end
    array_sort = sort(array,2, 'descend');
    maxs = array_sort(:,1:numb_highest);
    med_highest(:,1) = median(maxs,2);
    
else
    disp ('limited to 2D application')
end

end

