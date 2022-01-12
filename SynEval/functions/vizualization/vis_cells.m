function [I_stack]=vis_cells(image,BW,varargin)
%default
showI = 1;
color = 'red';

if nargin >= 3
    color = varargin{1};
end

if nargin >= 4
    showI = varargin{2};
end

if nargin == 5
    titel = varargin{3};
    figure ('Name',titel)
elseif showI == 0
else
    figure
end



maxi = max(image(:));
mini = min(image(:));
I_stack = (image-mini)/(maxi-mini);
BW_log = logical(BW);
I_stack(BW_log) = 2;
%I_stack = repmat(I_stack,[1 1 3]);
I_stack_temp = I_stack;

I_red = I_stack;
I_green = I_stack;
I_blue = I_stack;

switch color
    case {'red','color1'}
        I_red(I_stack_temp == 2) = 1;
        I_green(I_stack_temp == 2) = 0.2;
        I_blue(I_stack_temp == 2) = 0.2;
    case {'green','color2'}
        I_red(I_stack_temp == 2) = 0/255;
        I_green(I_stack_temp == 2) = 252/255;
        I_blue(I_stack_temp == 2) = 0;
    case {'orange','color3'}
        I_red(I_stack_temp == 2) = 255/255;
        I_green(I_stack_temp == 2) = 164/255;
        I_blue(I_stack_temp == 2) = 0;
    case 'blue'
        I_red(I_stack_temp == 2) = 0/255;
        I_green(I_stack_temp == 2) = 0/255;
        I_blue(I_stack_temp == 2) = 255/255;
    case 'magenta'
        I_red(I_stack_temp == 2) = 255/255;
        I_green(I_stack_temp == 2) = 0/255;
        I_blue(I_stack_temp == 2) = 255/255;
    otherwise
        errormessage
end

I_stack = cat(3,I_red, I_green, I_blue);

if showI == 1
    imshow(I_stack);
end
end