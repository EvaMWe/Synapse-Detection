% Preprosessing - version 1.0
% first version to load to github
% summarize PreProcessing (former version restricted to staining analysis)
% and PreProcessing_Imaging_ (including discard cell bodies)
%
%subsequent threshold calcualtion and S/N calcualtion are removed, since
%not necessary anymore due to iterative thresholding in subsequent
%celldetection process

function I = PreProcessing(img, discard, denoise, show_control_image)
%backgroundsubtraction
%normalization


img=double(img); 
if show_control_image == 1
figure(1)
imshow(img,'DisplayRange',[min(img(:)) max(img(:))])
end

%% denoise
if denoise == 1
    h = fspecial('gaussian',3, 1); %establish gaussian filter
    img_filtered = imfilter(img,h);
else
    img_filtered = img;
end

%% background Image
K_n = createKernel(21);
blurred = img_filtered;
blurred = conv2(blurred,K_n,'same');
template= blurred;

%% DISCARD CELLBODIES:
if discard == 1
SE = strel('disk',10,0);
Values = img_filtered(:);
[binData,intensity] = histcounts(Values);
%from the left side
binDataInv = binData(end:-1:1);
L = 1:length(binData);
CH = cumsum(binDataInv);
CMA = CH./L; %+1
CMAdiff = diff(CMA);        %+2
CMAlog= CMAdiff >= - 0.5;       %like the slope, looking for plateau
idx = find(CMAlog); 
Sbright = idx(1);
Sbright = length(binData) - Sbright; %add 3 due to starting from 2nd value and diff and take the end edge of bar 


%from the dark side
start = idx(end);
trace = CMA(start:-1:Sbright);
tracediff = diff(trace);
tracelog = tracediff >= - 0.5;
stop_ = find(tracelog,1);
if isempty(stop_)
    stop_ = length(trace);
end
Sdark = (length(binData) - start)+ stop_;

S = round((Sdark + Sbright)/2);
IntS = intensity(S);

maskImg = false(size(img));
maskImg(img_filtered>IntS) = 1;
maskDil = imdilate(maskImg,SE);



img_sub = img_filtered;
temp = template;
img_sub(maskDil) = 0;
temp(~maskDil) =0;
ImageSub = img_sub + temp;
% background_value = BackSubst_Gauss(img_filtered,0);
% [imageSub,BW] = subractCellBodies(img_filtered,background_value);
%apply kernel

% BG = ImageSub;
% BG = conv2(BG,K_n,'same');
else
    ImageSub = img_filtered;

end
I = ImageSub - template;
I(I<0) = 0;

end






function K_n = createKernel(siz)
%% create kernel
sequence = abs(-siz+1:-1);
seq2 = 2:siz-1;
sequence = [sequence, seq2];
kernel1 = repmat(sequence,length(sequence),1);
kernel2 = repmat(sequence',1,length(sequence));
K = sqrt(kernel1.*kernel1 + kernel2.*kernel2);
kn = sum(sum(K));
K_n = K./kn;
end
