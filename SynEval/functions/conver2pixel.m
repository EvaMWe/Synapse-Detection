function [MinArea,MaxArea] = convert2pixel(Settings) %unit is µm


%default
PixelSizeCamera = Settings.CamPixelSize;
binning = Settings.Binning;
magnification = Settings.Magnification;
MinimumSize = Settings.MinSize;
MaximumSize = Settings.MaxSize; 


%% calcualtions
sizePixel = (PixelSizeCamera/magnification)*binning;
% diameter in pixel
MinPixel = ceil(MinimumSize/sizePixel);
MaxPixel = ceil(MaximumSize/sizePixel);

%area in pixel
MinArea = floor((MinPixel/2)^2*pi);
MaxArea = floor((MaxPixel/2)^2*pi);
end
