%% this is a function, that takes 2 puncta from two differente images, crop them to similar size, for 
%calculation of correlation coefficients

function [puncta1, puncta2] = createObjectArray(PixelList1, PixelList2, PixelIdxList1,PixelIdxList2,image1,image2)
% PixelList: x y coordinates
xmin = min([PixelList1(:,1) PixelList2(:,1)]);
xmax = max([PixelList1(:,1) PixelList2(:,1)]);
ymin = min([PixelList1(:,2) PixelList2(:,2)]);
ymax = max([PixelList1(:,2) PixelList2(:,2)]);

section = zeros(ymax-ymin+1,xmax-xmin+1);
%transfer of marker 1
coordTemp1 = PixelList1;
coordTemp1(:,1) = coordTemp1(:,1) - xmin+1;
coordTemp1(:,2) = coordTemp1(:,2) - ymin+1;
sz = size(section);
linTemp1 = sub2ind(sz,coordTemp1(:,2),coordTemp1(:,1)); %!! attention: x is column, y is row
values = image1(PixelIdxList1);
puncta1 = section;
puncta1(linTemp1)=values;

%transfer of marker 1
coordTemp2 = PixelList2;
coordTemp2(:,1) = coordTemp2(:,1) - xmin+1;
coordTemp2(:,2) = coordTemp2(:,2) - ymin+1;
sz = size(section);
linTemp2 = sub2ind(sz,coordTemp2(:,2),coordTemp2(:,1)); %!! attention: x is column, y is row
values = image1(PixelIdxList2);
puncta2 = section;
puncta2(linTemp2)=values;

end