function ImgCol = convGrey(image, varargin)

%default
col = 'red';

if nargin == 2
    col = varargin{1};
end

maxi = max(image(:));
mini = min(image(:));

%[m n r] = size(image);
%just for visulaization purpose
Inorm= (image-mini)/(maxi-mini);
colImg = cat(3, Inorm, Inorm ,Inorm);

switch col
    case 'red'
        ImgCol = colImg;
        ImgCol(:,:,3) = 0;
        ImgCol(:,:,2) = 0;
    case 'green'
        ImgCol = colImg;
        ImgCol(:,:,1) = 0;
        ImgCol(:,:,3) = 0;
    otherwise
        errormessage
end

figure
imshow(ImgCol);
end