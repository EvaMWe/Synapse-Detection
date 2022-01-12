%% this is a function, that takes 2 puncta from two differente images, crop them to similar size,
% calculation of correlation
% 1. get range of coordinates (min and max) to span a template/section
% 2. take coordinates of puncta and transfer them to new coordinate system.
%    (x-xmin+1)
% 3. get values by converting subscipts o linear indices and fill them with the values derived from image and
%    former linear indices

function [ Overlap] = getOverlap(PixelList1, PixelList2, PixelIdxList1,PixelIdxList2,image)
% PixelList: x y coordinates
xmin = min([PixelList1(:,1); PixelList2(:,1)]);
xmax = max([PixelList1(:,1); PixelList2(:,1)]);
ymin = min([PixelList1(:,2); PixelList2(:,2)]);
ymax = max([PixelList1(:,2); PixelList2(:,2)]);

section = zeros(ymax-ymin+1,xmax-xmin+1);
%transfer of marker 1
coordTemp1 = PixelList1;
coordTemp1(:,1) = coordTemp1(:,1) - xmin+1;
coordTemp1(:,2) = coordTemp1(:,2) - ymin+1;
sz = size(section);
linTemp1 = sub2ind(sz,coordTemp1(:,2),coordTemp1(:,1)); %!! attention: x is column, y is row
values1 = image(PixelIdxList1);
puncta1 = section;
puncta1(linTemp1)=values1;

%transfer of marker 1
coordTemp2 = PixelList2;
coordTemp2(:,1) = coordTemp2(:,1) - xmin+1;
coordTemp2(:,2) = coordTemp2(:,2) - ymin+1;
sz = size(section);
linTemp2 = sub2ind(sz,coordTemp2(:,2),coordTemp2(:,1)); %!! attention: x is column, y is row
values2 = image(PixelIdxList2);
puncta2 = section;
puncta2(linTemp2)=values2;

%% overlap
logic1 = puncta1 > 0;
logic2 = puncta2 > 0;
sumLog = logic1+logic2;
Overlap = sum(sum(sumLog == 2));


end