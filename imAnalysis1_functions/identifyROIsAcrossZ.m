function [CaROImasks,userInput,ROIorders] = identifyROIsAcrossZ(reg_Stacks,userInput,UIr,numZplanes)

stackAVs = cell(1,length(reg_Stacks));
CaROImasks = cell(1,length(reg_Stacks));
ROIorders = cell(1,length(reg_Stacks));
for Z = 1:length(reg_Stacks)               
    %create average images per Z plane 
    stackAVs{Z} = mean(reg_Stacks{Z},3); 
    %scale the images to be between 0 to 1 
    stackAVs{Z} = stackAVs{Z} ./ max(stackAVs{Z}(:));

    CAroiGen = 1;
    while CAroiGen == 1
        [userInput,nm1BW2,UIr,CAroiGen] = createCaROIs(userInput,stackAVs{Z},UIr,Z);           
    end 
    CaROImasks{Z} = nm1BW2; 

    %figure out the order that terminal ROIs are looked at to match ROI
    %with data 
    ROIorders{Z} = bwlabel(nm1BW2);
end

meanPixIntArray = cell(1,size(reg_Stacks,2));
%apply mask to data and get the average pixel intensity of each terminal ROI 
for zStack = 1:size(reg_Stacks,2)
    for frame = 1:length(reg_Stacks{zStack})
        stats = regionprops(CaROImasks{zStack},reg_Stacks{zStack}(:,:,frame),'MeanIntensity');
            for ROI = 1:size(stats,1)
                meanPixIntArray{zStack}(ROI,frame) = stats(ROI).MeanIntensity;
            end 
    end   
end 

%compare the GCaMP signal of each Ca ROI with eachother 
%corr2 gives you a single correlation value for the entire vectors
%being compared 
ROIcorArray = cell(1,size(reg_Stacks,2));
for zStack = 1:size(reg_Stacks,2)
    for ROI = 1:size(meanPixIntArray{zStack},1) 
         for ROIcor = 1:size(meanPixIntArray{zStack},1)
            ROIcorArray{zStack}(ROI,ROIcor) = corr2(meanPixIntArray{zStack}(ROI,:),meanPixIntArray{zStack}(ROIcor,:));
         end
    end
end 

%this updates the Ca ROIs based on their activity patterns and
%location (2D) 
count = 1;
R = cell(1,count);
C = cell(1,count);
for Z = 1:size(reg_Stacks,2)
    for row = 1:size(meanPixIntArray{Z},1) 
        for col = 1:size(meanPixIntArray{Z},1) 
            if ROIcorArray{Z}(row,col) >= 0.75 && ROIcorArray{Z}(row,col) ~= 1 
                %find correlated ROIs 
                R{count} = ROIorders{Z} == row;
                C{count} = ROIorders{Z} == col;

                [CaROImasks,count] = combineCaROIs(CaROImasks,count,Z,R,C);                    
            end
        end 
    end 
end

%reapply mask to data and get the average pixel intensity of each terminal ROI 
clear meanPixIntArray 
for zStack = 1:size(reg_Stacks,2)
    for frame = 1:length(reg_Stacks{zStack})
        stats = regionprops(CaROImasks{zStack},reg_Stacks{zStack}(:,:,frame),'MeanIntensity');
            for ROI = 1:size(stats,1)
                meanPixIntArray{zStack}(ROI,frame) = stats(ROI).MeanIntensity;
            end 
    end   
end        

%3D caROI comparisons 
ROIcorArrayAcrossZ = cell(1,numZplanes);
for Z = 1:numZplanes
    for caROI = 1:size(meanPixIntArray{Z},1) 
        if Z < numZplanes
            %FIRST STEP = PICK ONE CAROI IN Z1 AND DO CROSS
            %CORRELATION OF THIS ROI WITH ALL ROIS IN Z2 

            for caROIZ2 = 1:size(meanPixIntArray{Z+1},1) 
                %compare the GCaMP signal of CaROI in Z with all
                %caROIs in Z+1
                %corr2 gives you a single correlation value for the entire vectors
                %being compared 
                ROIcorArrayAcrossZ{Z}(caROI,caROIZ2) = corr2(meanPixIntArray{Z}(caROI,:),meanPixIntArray{Z+1}(caROIZ2,:));
            end 
        end 
    end 
end 

%find centroid of each ROI
centroids = cell(1,numZplanes);
for Z = 1:numZplanes
    for caROI = 1:size(meanPixIntArray{Z},1) 
        centroids{Z} = regionprops(CaROImasks{Z},'centroid');   
    end 
end 

corrIDXs = cell(1,size(ROIcorArrayAcrossZ,2));
%find indices of correlated ROIs across Z 
for z = 1:size(ROIcorArrayAcrossZ,2)
     [x, y] = find(ROIcorArrayAcrossZ{z} > 0.01);
     if isempty(x) == 0 
         corrIDXs{z}(:,1) = x;
     end 
     if isempty(y) == 0
         corrIDXs{z}(:,2) = y;
     end 
end 

centroidLoc = cell(1,size(corrIDXs,2));
centroidLoc2 = cell(1,size(corrIDXs,2));
CorrROIdists = cell(1,size(corrIDXs,2));
for z = 1:size(corrIDXs,2) 
    for corrROI = 1:size(corrIDXs{z},1)
        %get the locations of correlated ROIs 
        centroidLoc{z}{corrROI} = centroids{z}(corrIDXs{z}(corrROI,1)).Centroid;
        centroidX = centroidLoc{z}{corrROI}(1,1);
        centroidY = centroidLoc{z}{corrROI}(1,2);

        centroidLoc2{z}{corrROI} = centroids{z+1}(corrIDXs{z}(corrROI,2)).Centroid;
        centroidX2 = centroidLoc2{z}{corrROI}(1,1);
        centroidY2 = centroidLoc2{z}{corrROI}(1,2);

       %compute distances between all correlated ROIs across Z
       CorrROIdists{z}(corrROI) = sqrt(((centroidX-centroidX2).^2) + ((centroidY-centroidY2).^2));               
    end 
    CorrROIdists{z} = CorrROIdists{z}';
    if isempty(CorrROIdists{z}) == 0
        corrIDXs{z}(:,3) = CorrROIdists{z};
    end 
end 

%finds ROIs that are correlated with more than one other ROI (the
%repeats) 
corrIDXnoRepeatsCol1 = cell(1,size(corrIDXs,2));
corrIDXnoRepeatsCol2 = cell(1,size(corrIDXs,2));
for z = 1:size(corrIDXs,2) 
    if isempty(corrIDXs{z}) == 0
        [noRepeatsCol1,~,~] = unique(corrIDXs{z}(:,1));
        [noRepeatsCol2,~,~] = unique(corrIDXs{z}(:,2));
        corrIDXnoRepeatsCol1{z} = noRepeatsCol1;
        corrIDXnoRepeatsCol2{z} = noRepeatsCol2;
    end 
end 


%remove repeating ROIs based on distances 
for z = 1:size(corrIDXs,2) 
    for uniqueVals = 1:length(corrIDXnoRepeatsCol1{z})
        %if an ROI is repeating col 1 
        if length(find(corrIDXs{z}(:,1) == corrIDXnoRepeatsCol1{z}(uniqueVals))) > 1 
            %find repeating ROI indices 
            repIDXs = find(corrIDXs{z}(:,1) == corrIDXnoRepeatsCol1{z}(uniqueVals));
            %get the distances 
            dists = zeros(1,length(repIDXs));
            for repIDs = 1:length(repIDXs)
                dists(repIDs) = corrIDXs{z}(repIDXs(repIDs),3);
            end 
            %get min distance 
            minDist = min(dists);

            %remove ROI Z-Z2 pairs that are not the min dist - this
            %blanks the data with rows of 0s 
            for repIDs = 1:length(repIDXs)
                if corrIDXs{z}(repIDXs(repIDs),3) ~= minDist
                    corrIDXs{z}(repIDXs(repIDs),:) = [0,0,0];
                end 
            end 
        end 
    end 
    
    for uniqueVals = 1:length(corrIDXnoRepeatsCol2{z})
        %if an ROI is repeating col 2 
        if length(find(corrIDXs{z}(:,1) == corrIDXnoRepeatsCol2{z}(uniqueVals))) > 1 
            %find repeating ROI indices 
            repIDXs = find(corrIDXs{z}(:,1) == corrIDXnoRepeatsCol2{z}(uniqueVals));
            %get the distances 
            for repIDs = 1:length(repIDXs)
                dists(repIDs) = corrIDXs{z}(repIDXs(repIDs),3);
            end 
            %get min distance 
            minDist = min(dists);

            %remove ROI Z-Z2 pairs that are not the min dist - this
            %blanks the data with rows of 0s 
            for repIDs = 1:length(repIDXs)
                if corrIDXs{z}(repIDXs(repIDs),3) ~= minDist
                    corrIDXs{z}(repIDXs(repIDs),:) = [0,0,0];
                end 
            end 
        end                
    end 
end 

 %below removes correlated ROIs across Z that are still too far
 %away to belong to the same cell 
 for z = 1:size(corrIDXs,2)
     for ROIpair = 1:size(corrIDXs{z},1)
         if corrIDXs{z}(ROIpair,3) > 2 %2 here is setting a distance threshold 
             corrIDXs{z}(ROIpair,:) = [0,0,0];
         end 
     end 
 end 

 %remove rows that are full of zeros 
 corrIDXs2 = cell(1,size(corrIDXs,2));
 for z = 1:size(corrIDXs,2)
     %this matrix is gives us the ROI pairs across Z!!! 
     corrIDXs2{z} = corrIDXs{z}(any(corrIDXs{z},2),:);
 end 

%relables ROIs that belong to the same terminal across Z        
for x = 1:length(CaROImasks)
    CaROImasks{x} = double(CaROImasks{x});
end 
count = 2;
for z = 1:size(corrIDXs2,2)
    if z == 1 
        for ROIpair = 1:size(corrIDXs2{z},1)                 
            CaROImasks{z}(ROIorders{z} == corrIDXs2{z}(ROIpair,1)) = count; 
            CaROImasks{z+1}(ROIorders{z+1} == corrIDXs2{z}(ROIpair,2)) = count;
            count = count+1; 
        end 
    elseif z > 1 
        for ROIpair = 1:size(corrIDXs2{z},1) 
            %if a Z2 ROI ind in Z2/3 pair is in a Z1/2 pair use that
            %ROIs count number in Z3 THIS MUST WORK FOR HOWEVER
            %MANY Z PLANES THERE ARE 
            if ismember(corrIDXs2{z}(ROIpair,1),corrIDXs2{z-1}(:,2)) == 1 
                CaROImasks{z+1}(ROIorders{z+1} == corrIDXs2{z}(ROIpair,2)) = unique(CaROImasks{z}(ROIorders{z} == corrIDXs2{z}(ROIpair,1)));  
            end 
        end 
    end 
end          

%figure out where to start label value for ROIs that don't span
%multiple Z planes 
maxZ = zeros(1,length(meanPixIntArray));
for Z = 1:length(meanPixIntArray)
    maxZ(Z) = max(max(CaROImasks{Z}));
end 
maxLabel = max(maxZ);

%relable all ROIs in CaROImasks that do not appear in multiple Z
%planes w/non 1 numbers 
maxLabel = maxLabel + 1;
for Z = 1:length(meanPixIntArray)
    for ROIorder = 1:size(meanPixIntArray{Z},1)
        if unique(CaROImasks{Z}(ROIorders{Z} == ROIorder)) == 1 
            CaROImasks{Z}(ROIorders{Z} == ROIorder) = maxLabel;
            maxLabel = maxLabel + 1;
        end 
    end 
end 


imagesc(CaROImasks{1});grid on;figure;imagesc(CaROImasks{2});grid on;figure;imagesc(CaROImasks{3});grid on;
end 