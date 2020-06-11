function [CaROImasks,count] = combineCaROIs(CaROImasks,count,Z,R,C)



%determine average size of the correlated terminals
for ROI = 1:length(R)
    %determine number of white pixels per row 
    white_pixR = flipud(sum(R{ROI} == 1, 2));                     
    %find max ROI width
    Max_white_pixR = max(white_pixR);
end 
for ROI = 1:length(C)
    %determine number of white pixels per row 
    white_pixC = flipud(sum(C{ROI} == 1, 2));                     
    %find max ROI width
    Max_white_pixC = max(white_pixC);
end 
ROIwidths(1) = Max_white_pixR; ROIwidths(2) = Max_white_pixC;
avROIwidth = mean(ROIwidths,2);

%find location of corralted ROIs 
[xR,yR] = find(R{count} == 1);
[xC,yC] = find(C{count} == 1);
clear Rinds Cinds 
Rinds(:,1) = xR;Rinds(:,2) = yR;
Cinds(:,1) = xC;Cinds(:,2) = yC;
count = count + 1;

%find shortest distance between correlated ROIs 
ROIdistances = zeros(size(Rinds,1),size(Cinds,1));
for ROI1 = 1:size(Rinds,1)
    for ROI2 = 1:size(Cinds,1)
        ROIdistances(ROI1,ROI2) = sqrt((Rinds(ROI1,1)-Cinds(ROI2,1)).^2 + (Rinds(ROI1,2)-Cinds(ROI2,2)).^2);
    end 
end 
minMinDist = min(min(ROIdistances));

%set distance threshold
if minMinDist < avROIwidth
    %dilate relevant ROIs 
    radius = round(minMinDist);
    decomposition = 0;
    se = strel('disk', radius, decomposition);
    Rdilated = imdilate(R{1},se);
    Cdilated = imdilate(C{1},se);

    %get indices of the dilated ROI
    clear RdilInds CdilInds
    [xRdil,yRdil] = find(Rdilated == 1);
    RdilInds(:,1) = xRdil;RdilInds(:,2) = yRdil;
    [xCdil,yCdil] = find(Cdilated == 1);
    CdilInds(:,1) = xCdil;CdilInds(:,2) = yCdil;

    %find where the dilated ROI and correlated ROI
    %overlap 
    for RdilIndsRow = 1:size(RdilInds,1)                          
        for CindsRow = 1:size(Cinds,1)
            closestCcheck(CindsRow,:,RdilIndsRow) = RdilInds(RdilIndsRow,:) == Cinds(CindsRow,:);                             
        end 
    end 
    for CdilIndsRow = 1:size(CdilInds,1)                              
        for RindsRow = 1:size(Rinds,1)
            closestRcheck(RindsRow,:,CdilIndsRow) = CdilInds(CdilIndsRow,:) == Rinds(RindsRow,:);                             
        end 
    end 

    %determine closest pixels based on overlap
    overlap = [1 1];
    count2 = 1;
    for RdilIndsRow = 1:size(RdilInds,1)
        [~, closestCind] = ismember(overlap, closestCcheck(:,:,RdilIndsRow), 'rows');
        if closestCind > 0 
            closestCpixLoc(count2,:) = Cinds(closestCind,:);
        end   
        count2 = count2+1;
    end 
    count3 = 1;
    for CdilIndsRow = 1:size(CdilInds,1)
        [~, closestRind] = ismember(overlap, closestRcheck(:,:,CdilIndsRow), 'rows');
        if closestRind > 0 
            closestRpixLoc(count3,:) = Rinds(closestRind,:);
        end   
        count3 = count3+1;
    end 
    
    %remove zeros from closestCpixLoc and closestRpixLoc 
    closestRpixLoc2 = closestRpixLoc(any(closestRpixLoc,2),:);
    closestCpixLoc2 = closestCpixLoc(any(closestCpixLoc,2),:);
    
    %remove repeating rows 
    closestRpixLoc3 = unique(closestRpixLoc2,'rows');
    closestCpixLoc3 = unique(closestCpixLoc2,'rows');
    
    %identify the in b/w pixels that will now be a part
    %of the mask/ROI 
    if size(closestRpixLoc3,1) > size(closestCpixLoc3,1)
        numPix = size(closestCpixLoc3,1);
    elseif size(closestRpixLoc3,1) < size(closestCpixLoc3,1)
        numPix = size(closestRpixLoc3,1);
    elseif size(closestRpixLoc3,1) == size(closestCpixLoc3,1)
        numPix = size(closestRpixLoc3,1);
    end 

    bwPixLoc = zeros(numPix,1);
    for bwPix = 1:numPix
        %x coordinate
        if closestRpixLoc3(bwPix,1) == closestCpixLoc3(bwPix,1)
            bwPixLoc(bwPix,1) = closestRpixLoc3(bwPix,1);
        elseif closestRpixLoc3(bwPix,1) ~= closestCpixLoc3(bwPix,1)
            bwPixLoc(bwPix,1) = (closestRpixLoc3(bwPix,1)+closestCpixLoc3(bwPix,1))/2;
        end 
        %y coordinate 
        if closestRpixLoc3(bwPix,2) == closestCpixLoc3(bwPix,2)
            bwPixLoc(bwPix,2) = closestRpixLoc3(bwPix,2);
        elseif closestRpixLoc3(bwPix,2) ~= closestCpixLoc3(bwPix,2)
            bwPixLoc(bwPix,2) = (closestRpixLoc3(bwPix,2)+closestCpixLoc3(bwPix,2))/2;
        end 
    end 
    
    %remove indices that are not unique integers 
    bwPixLoc = floor(bwPixLoc); bwPixLoc = unique(bwPixLoc,'rows');

    %change the original mask to update the ROIs 
    for bwPix = 1:size(bwPixLoc,1)
        CaROImasks{Z}(floor(bwPixLoc(bwPix,1)),floor(bwPixLoc(bwPix,2))) = 1; 
    end 
                
end 