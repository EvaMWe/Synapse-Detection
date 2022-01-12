%% function to transfer coodinates
function idxList = coordTrans(idxListtrans,sizSec,sizIm,xmin,ymin)
idxListID = idxListtrans(:,1);
[r_t,c_t] = ind2sub(sizSec,idxListID);
x = c_t+xmin-3;
y = r_t+ymin-3;

neg_x = find(x<=0);
neg_y = find(y<=0);
neg = [neg_x;neg_y];
x(neg) = [];
y(neg) = [];
idxList = sub2ind(sizIm,y,x);
valList = idxListtrans(:,2);
valList(neg) = [];
idxList = [idxList valList];
end