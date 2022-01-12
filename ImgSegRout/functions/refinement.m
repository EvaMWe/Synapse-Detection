function parameter = refinement (parameter)     

X = parameter.centerX;
Y = parameter.centerY;
img = parameter.img;
area = parameter.area;
m0 = parameter.m0;
idx = parameter.idx;
idy = parameter.idy;
w = parameter.radius;
mask = parameter.mask;

errX = sum(sum(area.*idx))/m0;
count = 0;
while abs(errX) > 0.5 && count <= 3
    count = count +1;
    X = X - sign(errX);
    if X-w <= 0 || X+w > size(img,2) || Y-w <= 0 || Y+w > size(img,1)
        break
    end
    area = img(Y-w:Y+w,X-w:X+w);    
    area(mask == 0) = 0;
    m0 = sum(sum(area));
    errX = sum(sum(area.*idx))/m0;    
end
    
errY = sum(sum(area.*idy))/m0;
count = 0;
while abs(errY) > 0.5 && count <= 3
    count = count +1;
    Y = Y - sign(errY);
    if Y-w <= 0 || Y+w > size(img,1) || X-w <= 0 || X+w > size(img,2)
        break
    end
    area = img(Y-w:Y+w,X-w:X+w);    
    area(mask == 0) = 0;
    m0 = sum(sum(area));
    errY = sum(sum(area.*idx))/m0;
end

parameter.area = area;
parameter.m0 = m0;
parameter.centerX = X;
parameter.centerY = Y;
        
end