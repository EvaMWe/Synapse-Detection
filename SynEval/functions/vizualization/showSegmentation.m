%Visualisieren der Daten

function I = showSegmentation (image_norm ,BW, regList)

SE = strel('disk',1,0);

circleMarks = false(size(BW));
nReg = length(regList);
 for r = 1:nReg
     template = false(size(BW));
     idx = regList(r).PixelIdxList;
     template(idx) = 1;
     sub = imerode(template,SE);
     circle = logical(template - sub);
     circleMarks(circle) = 1;
 end
 
 I = vis_cells(image_norm,circleMarks, 'red',0);
end

%Normalizing
