function [SN,signal,noise] =calc_sig2noise(image)
ImgVec = image(:);
ImgVec = double(ImgVec);
[counts,bins] = hist(ImgVec,500);
%hist(ImgVec,50) falls mans doch mal sehen will

% get area for noise calcualtion
k = 0;
pth = pth/100;
while sum(counts(end-k:end))< pth*sum(counts) %threshold finden
    k = k+1;
end
level_noise=bins(end-k+1);
BG = ImgVec;
BG(BG> level_noise) = 0;
noise = std(BG);

% get area for signal  calculation
pth = 0.9;
while sum(counts(end-k:end))< pth*sum(counts) %threshold finden
    k = k+1;
end
level_sig=bins(end-k+1);
Sig = ImgVec;
Sig(Sig < level_sig) = 0;
signal = mean(Sig);

SN = signal/noise;
end