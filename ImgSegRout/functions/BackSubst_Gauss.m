function [mean_img, level] = BackSubst_Gauss(img, show_control_image)
%Prozedur zum Hintergrundeliminieren für die DFG-Messungen

pth = 0.1;
img=double(img); 

if show_control_image == 1
figure(1)
imshow(img,'DisplayRange',[min(img(:)) max(img(:))])
end

%artificial blurring
h = fspecial('gaussian',5, 1); %establish gaussian filter
img_filtered = imfilter(img,h); 

if show_control_image == 1
    figure(2)
    imshow(img_filtered,'DisplayRange',[min(img(:)) max(img(:))])
    title('Image after gaussian filter')
end

ImgVec = img_filtered(:);
ImgVec = double(ImgVec);
[counts,bins] = hist(ImgVec,500);
%hist(ImgVec,50) falls mans doch mal sehen will
k = 0;
pth = pth/100;
while sum(counts(end-k:end))< pth*sum(counts) %threshold finden
    k = k+1;
end
level=bins(end-k+1);


%Background-Image als Grau-Bild
BG = img_filtered;
BG(BG > level) = 0;

if show_control_image == 1
figure(3)
imshow(BG,'DisplayRange',[min(img(:)) max(img(:))])
end

BG(BG==0) = NaN;
mean_img = nanmean(BG,1); 
mean_img = nanmean(mean_img);

end