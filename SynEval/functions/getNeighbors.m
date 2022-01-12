%This function returns the neighbours around one specific pixel; 
%Thereby w defines the area around the pixel; w = 1, just direct neighbours
% connected to the pixel;
%1. List containing the linear indices of the neighbors
%2. Values belonging to the indices

function [idx_list, values] = getNeighbors(section, centroids, w)

siz = size(section);
n_values = (w*2+1)^2;
n_rows = size(section,1);
n_col = size(section,2);
% get idx from 8 neighbors
map_rows =repmat((1:n_rows)',1,n_col);
map_col = repmat(1:n_col,n_rows,1);

cen_row = map_rows(centroids);
cen_col = map_col(centroids);

mask_row = repmat((-w:w)',1,2*w+1);
rows  = mask_row + cen_row;
mask_col = repmat((-w:w),2*w+1,1);
columns = mask_col + cen_col;

vec_row = reshape(rows,n_values,1);
vec_col = reshape(columns,n_values,1);

idx_list = sub2ind(siz,vec_row,vec_col);
idx_list(idx_list == centroids) = [];
values = section(idx_list);

end