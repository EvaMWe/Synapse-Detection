function [Preprocessing, Iteration, SizeSettings] = configuration()
Preprocessing.Discard = [];
Preprocessing.Denoise = [];
Preprocessing.ShowImage = [];

Iteration.it = [];

SizeSettings.CamPixelSize = [];
SizeSettings.Magnification = [];
SizeSettings.Binning = [];
SizeSettings.MaxSize = [];
SizeSettings.MinSize = [];

end
