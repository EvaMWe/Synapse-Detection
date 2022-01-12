%update: 06.01.2022
function  regList = precising_woList(ListOfAll,numbC, sizImg,  h, image)
counter = 1;
regList = zeros(400,numbC*5);
template = zeros(sizImg);
SE = strel('disk',1,0);
SE2 = strel('disk',2,0);
for n = 1:numbC
    area = ListOfAll((ListOfAll(:,2) == n),1);
    areaRest = ListOfAll((ListOfAll(:,2) ~= n &(ListOfAll(:,2) ~= numbC+1)),1);
    %% refine
    if length(area) <= 6
        continue
    end
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
    end
    %%
    
    regCleaned_ = cleanRg(area,sizImg,SE);
    template_sel = zeros(sizImg);
    template_sel(regCleaned_) = 1;
    
    %% check if there are more than one regions
    select = bwconncomp(template_sel,4);
    if select.NumObjects > 1
        numObj = select.NumObjects;
        
        for nO = 1:numObj           
            PixelIdxList = select.PixelIdxList{1,nO};            
            if length(PixelIdxList) > 20
                regCleaned = cleanRg(PixelIdxList,sizImg,SE);
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
            elseif length(PixelIdxList) < 6
                 continue
            else
                regCleaned = PixelIdxList;
                %to be sure
                templClean = zeros(sizImg);
                templClean(regCleaned) =1;
                sel = bwconncomp(templClean,4);
                if length(sel.PixelIdxList) > 1
                    siz = cellfun(@length,sel.PixelIdxList);
                    [~,maxIdx]= max(siz);
                    regCleaned = sel.PixelIdxList{1,maxIdx};
                end
            end
            
            if length(regCleaned) < 6
                continue
            end
            
            regList(1:length(regCleaned),counter) = regCleaned;           
            counter = counter +1;
        end
    else
        if length(area) > 20
            regCleaned = cleanRg(area,sizImg,SE);
        else
            regCleaned = area;
        end
         %to be sure
         templClean = zeros(sizImg);
         templClean(regCleaned) =1;
         sel = bwconncomp(templClean,4);
         if length(sel.PixelIdxList) > 1
             siz = cellfun(@length,sel.PixelIdxList);
             maxIdx = find(siz>6);
             nbIdx = length(maxIdx);
             for id = 1:nbIdx
                 regCleaned = sel.PixelIdxList{1,maxIdx(id)};
                 regList(1:length(regCleaned),counter) = regCleaned;
                 counter = counter +1;
             end
         else
             regList(1:length(regCleaned),counter) = regCleaned;
             counter = counter +1;
         end
         
    end
    
end
regList = regList(:,1:counter-1);
end

%clean regions from tails
function regCleaned = cleanRg(PixelIdxList,sizImg,SE)
BWcur = zeros(sizImg);
BWcur(PixelIdxList) = 1;
BWdisc = imerode(BWcur,SE);
BWdisc = imdilate(BWdisc,SE);
BWdisc = imfill(BWdisc,4,'holes');
regCleaned = find(BWdisc ==1);
end
