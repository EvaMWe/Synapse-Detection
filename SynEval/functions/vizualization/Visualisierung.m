%Visualisieren der Daten

function Visualisierung (image_norm ,BW)
vis_cells(image_norm ,BW);

SE = strel('disk',1,0);
BW_small = imerode(BW,SE);

transparent(filtered,BW_ws_cleared,BW_small)


circleMarks = false(size(BW_ref));
 for r = 1:37
     template = false(size(BW_ref));
     idx = reg_ref(r).PixelIdxList;
     template(idx) = 1;
     sub = imerode(template,SE);
     circle = logical(template - sub);
     circleMarks(circle) = 1;
 end
vis_cells(croppedI,circleMarks, 'red',1);

%Normalizing
