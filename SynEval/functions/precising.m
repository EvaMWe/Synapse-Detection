function[counter, regList] = precising(regList,ListOfAll,numbC, sizImg, counter, h, image, minArea)
template = zeros(sizImg);
SE = strel('disk',1,0);
SE2 = strel('disk',2,0);
for n = 1:numbC
    area = ListOfAll((ListOfAll(:,2) == n),1);
    if length(area) < minArea
        continue
    end
    areaRest = ListOfAll((ListOfAll(:,2) ~= n &(ListOfAll(:,2) ~= numbC+1)),1);
    %% refine: search for neighboring pixel to include to puncta
    if exist('h','var') && exist('image', 'var')
        template(area) = 1;
        templatedil = imdilate(template, SE2);
        templateEnd = templatedil-template; %thats the border
        %templateEnd(image(logical(templateEnd)) < h) = 0;
        templateEnd = templateEnd & image > h;
        above = find(templateEnd == 1);
        area = [above; area];
        area = area(~ismember(area,areaRest)); %discard pixel that are already in other regions of this blobb
        ListOfAll(ismember(ListOfAll(:,1),area),2) = n;
        clear 'above'        
    end
    if length(area) < minArea 
        continue
    end
    %%
    template_sel = zeros(sizImg);
    template_sel(area) = 1;
    template_sel = imfill(template_sel,4,'holes');
    %%check if there are more than one regions
    select = bwconncomp(template_sel,4);
    if select.NumObjects > 1
        numObj = select.NumObjects;
        for nO = 1:numObj
            PixelIdxList = select.PixelIdxList{1,nO};
            regCleaned = cleanRg(PixelIdxList,sizImg,SE);
            if length(regCleaned) < minArea
                continue
            end
            % check for separate regions again (that developed due to
            % cleanRag
            templClean = zeros(sizImg);
            templClean(regCleaned) =1;
            sel = bwconncomp(templClean,4);
            if length(sel.PixelIdxList) > 1
                siz = cellfun(@length,sel.PixelIdxList);
                [~,maxIdx]= max(siz);
                regCleaned = sel.PixelIdxList{1,maxIdx};
            end
            regList(counter).PixelIdxList = regCleaned;
            regList(counter).Area = length(regCleaned);
            [r,c] = ind2sub(sizImg,regCleaned);
            regList(counter).PixelList = [r c];
            counter = counter+1;
        end
    else
        regCleaned = cleanRg(area,sizImg,SE);
        % check for separate regions again (that developed due to
        % cleanRag
        templClean = zeros(sizImg);
        templClean(regCleaned) =1;
        sel = bwconncomp(templClean,4);
        if length(sel.PixelIdxList) > 1
            siz = cellfun(@length,sel.PixelIdxList);
            [~,maxIdx]= max(siz);
            regCleaned = sel.PixelIdxList{1,maxIdx};
        end
        regList(counter).PixelIdxList = regCleaned;
        regList(counter).Area = length(regCleaned);
        [r,c] = ind2sub(sizImg,regCleaned);
        regList(counter).PixelList = [r c];
        counter = counter+1;
    end
    
end
end

%clean regions from tails
function regCleaned = cleanRg(PixelIdxList,sizImg,SE)
        BWcur = zeros(sizImg);
        BWcur(PixelIdxList) = 1;
        BWdisc = imerode(BWcur,SE);
        BWdisc = imdilate(BWdisc,SE);
        regCleaned = find(BWdisc ==1);
end
        