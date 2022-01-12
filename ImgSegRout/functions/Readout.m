%
%input: image stack: containing the frames
%       region properties: containing the coordinates for region of
%       interest
%output: restored traces from individual region of interest (normalized and
%        background subtracted)
%        Normalization to the unstimulated fluorescence value at the
%        beginning of the experiment

function [data, backgroundTrace, regionNb] = Readout(pixelIdx, regionNb, stack)

%calculation of the background trace


frameNb=size(stack,3);
backgroundTrace = zeros(1,frameNb);
for frame = 1:frameNb
    img = stack(:,:,frame);
    [~,background_value] = BackgroundSubtraction_Roling(img,21);
    backgroundTrace(frame) = background_value;
end
backgroundTrace = smooth(backgroundTrace)';

%[BG_matrix] = BackSubst_gradient(stack,pixelIdx);
%BG_matrix_ = arrayfun(@(M(i,:)) smooth(BG_matrix(i)),1:792);

%readout data
data = zeros(regionNb,frameNb);
for frame = 1:frameNb
    img = stack(:,:,frame);
    for region = 1:regionNb
        pxl = pixelIdx(region,1).PixelIdxList;
        if pxl == 1
            continue
        end
        data(region,frame) = mean(img(pxl));
    end
end

data = data(any(data,2),:);
regionNb = size(data,1);


%disp 'read-out finished'
end