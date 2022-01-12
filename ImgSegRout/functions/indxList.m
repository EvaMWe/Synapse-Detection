%get index list
% input:
% siz = image size, necessary for computing the indices out of subscripts
% centerC = subscipt of column
% centerR = subscript of row
function [pxlIndList] = indxList(siz, centerC, centerR, radius)


idx = -radius:1:radius;     % index vector
d = 2*radius+1;         % diameter
mask1 = repmat(idx',1,d);
mask2 = repmat(idx,d,1);
mask = zeros(d,d);
maskTemp = mask1.^2+mask2.^2;

mask(maskTemp <= radius*radius) = 1;

idC = centerC-radius:centerC+radius;
idR = centerR-radius:centerR+radius;

if sum(idC <= 0) ~= 0 || sum(idR <= 0) ~= 0
    pxlIndList = 0;
    return
end

idC = repmat(idC,2*radius+1,1);
idC = idC(:)';

idR = repmat(idR,1,2*radius+1);

if any(min(idR(:)) < 1) || any(max(idR(:)) > siz(1)) || any(min(idC(:)) < 1) || any(max(idC(:)) > siz(2))
   pxlIndList = 0;
  
else
    pxlIndList = sub2ind(siz,idR,idC);
    mat = vec2mat(pxlIndList,d);
    mat(mask == 0) = 0;
    pxlIndList = mat(:);
    pxlIndList(pxlIndList==0) = [];
end



end