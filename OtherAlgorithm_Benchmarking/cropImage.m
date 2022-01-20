%cropped I: cropped + BG subtracte, noise reduced, normalized
%IBG: cropped + BG subtracted, noise reduced

function [croppedI] = cropImage(x1,x2,y1,y2,varargin)
if nargin >= 5
    I = varargin{1};
else
    [imageFile, imagePath] = uigetfile('*.*','select image to crop');
    imageName = fullfile(imagePath,imageFile);
    if ~iscell(imageName)
        imageName = cellstr(imageName);
    end
    I = double(imread(imageName{1}));
end
    [IBG] = PreProcessing(I, 0,1,0);
    
    I_norm = (IBG-min(IBG(:)))/ max(IBG(:))-min(IBG(:));

croppedI = I_norm(y1:y2,x1:x2);

%croppedI = I_norm(y1:y2,x1:x2);
imshow(croppedI)
end