function [Preprocessing, Iteration, SizeSettings, FrameRange, Func] = configuration_ImgRegSeg()
Preprocessing.Discard = [];
Preprocessing.Denoise = [];
Preprocessing.ShowImage = [];

Iteration.it = [];

SizeSettings.CamPixelSize = [];
SizeSettings.Magnification = [];
SizeSettings.Binning = [];
SizeSettings.MaxSize = [];
SizeSettings.MinSize = [];

FrameRange.Frame1 = [];
FrameRange.Frame2 = [];

Func.Align = [];
Func.bleaching = [];

end
