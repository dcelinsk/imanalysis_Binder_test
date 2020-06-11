% %% register images 
clear AVsortedData AVwheelData indices 
% [regStacks,userInput,~,state_start_f,state_end_f,vel_wheel_data,TrialTypes,HDFchart] = imRegistration2(userInput);

%% set what data you want to plot 
[dataParseType] = getUserInput(userInput,'How many seconds before the stimulus starts do you want to plot?');
if dataParseType == 0 
    [sec_before_stim_start] = getUserInput(userInput,'What data do you need? Peristimulus epoch = 0. Stimulus epoch = 1.');
    [sec_after_stim_end] = getUserInput(userInput,'How many seconds after stimulus end do you want to plot?');
end 

%% set up analysis pipeline 
[VsegQ] = getUserInput(userInput,'Do you need to measure vessel width? Yes = 1. No = 0.');
[pixIntQ] = getUserInput(userInput,'Do you need to measure changes in pixel intensity? Yes = 1. No = 0.');
if pixIntQ == 1
    [CaQ] = getUserInput(userInput,'Do you need to measure changes in calcium dynamics? Yes = 1. No = 0.');
    [BBBQ] = getUserInput(userInput,'Do you need to measure changes in BBB permeability? Yes = 1. No = 0.');
end 
[cumStacksQ] = getUserInput(userInput,'Do you want to generate cumulative pixel intensity stacks? Yes = 1. No = 0.');

%% select registration method that's most appropriate for making the dff and cum pix int stacks 
[volIm] = getUserInput(userInput,'Is this volume imaging data? Yes = 1. Not = 0.');
if cumStacksQ == 1 || pixIntQ == 1 
    if volIm == 0
        regTypeDim = 0; 
    elseif volIm == 1 
        [regTypeDim] = getUserInput(userInput,'What registration dimension is best for pixel intensity analysis? 2D = 0. 3D = 1.');
    end 
    [regTypeTemp] = getUserInput(userInput,'What registration template is best for pixel intensity analysis? red = 0. green = 1.');
    [reg__Stacks] = pickRegStack(regStacks,regTypeDim,regTypeTemp);
    [reg_Stacks] = backgroundSubtraction2(reg__Stacks,BG_ROIboundData);
end
%% select registration method that's most appropriate for vessel segmentation 
if VsegQ == 1
    [regTypeDimVesSeg] = getUserInput(userInput,'What registration dimension is best for vessel segmentation? 2D = 0. 3D = 1.');
    [regTypeTempVesSeg] = getUserInput(userInput,'What registration template is best for vessel segmentation? red = 0. green = 1.');    
    [reg__StacksVesSeg] = pickRegStack(regStacks,regTypeDimVesSeg,regTypeTempVesSeg);
    [reg_Stacks] = backgroundSubtraction2(reg__StacksVesSeg,BG_ROIboundData);
end 

%% make cumulative, diff-cumulative, and DF/F stacks to output for calcium and BBB perm analysis 
[FPS] = getUserInput(userInput,"FPS"); 
[numZplanes] = getUserInput(userInput,"How many Z planes are there?");
% if cumStacksQ == 1  
%     [dffDataFirst20s2,CumDffDataFirst20s2,CumData2] = makeCumPixWholeEXPstacks(FPS,reg_Stacks,numZplanes,sec_before_stim_start);
% end

%% make sure state start and end frames line up 
[state_start_f,state_end_f,TrialTypes] = makeSureStartEndTrialTypesLineUp(reg_Stacks,state_start_f,state_end_f,TrialTypes,numZplanes);

%% resample velocity data by trial type 
% Dmitrijs
% if length(vel_wheel_data)*size(reg_Stacks{1},3) > 2^31
%     ResampedVel_wheel_data1 = resample(vel_wheel_data,(length(vel_wheel_data)/2),length(vel_wheel_data));
%     ResampedVel_wheel_data = resample(ResampedVel_wheel_data1,size(reg_Stacks{1},3),length(ResampedVel_wheel_data1));
% elseif length(vel_wheel_data)*size(reg_Stacks{1},3) < 2^31
%     ResampedVel_wheel_data = resample(vel_wheel_data,size(reg_Stacks{1},3),length(vel_wheel_data));
% end 

%% get rid of frames/trials where registration gets wonky 
%EVENTUALLY MAKE THIS AUTOMATIC INSTEAD OF HAVING TO INPUT WHAT FRAME THE
%REGISTRATION GETS WONKY 
UIr = size(userInput,1)+1;
cutOffFrameQ = input('Does the registration ever get wonky? Yes = 1. No = 0. ');  userInput(UIr,1) = ("Does the registration ever get wonky? Yes = 1. No = 0."); userInput(UIr,2) = (cutOffFrameQ); UIr = UIr+1;

if cutOffFrameQ == 1 
    cutOffFrame = input('Beyond what frame is the registration wonky? ');  userInput(UIr,1) = ("Beyond what frame is the registration wonky?"); userInput(UIr,2) = (cutOffFrame); UIr = UIr+1;
    if pixIntQ == 1 
        reg___Stacks = reg_Stacks; clear reg_Stacks; 
        reg_Stacks = cell(1,numZplanes);
        for zStack = 1:numZplanes
            reg_Stacks{zStack} = reg___Stacks{zStack}(:,:,1:cutOffFrame);
        end 
    end 
    
    if VsegQ == 1
        reg___StacksVesSeg = reg_Stacks; clear reg_Stacks; 
        reg_Stacks = cell(1,numZplanes);
        for zStack = 1:numZplanes
            reg_Stacks{zStack} = reg___StacksVesSeg{zStack}(:,:,1:cutOffFrame);
        end 
    end 
    
    if cumStacksQ == 1 
        reg___Stacks = reg_Stacks; clear reg_Stacks; 
        reg_Stacks = cell(1,numZplanes);
        for zStack = 1:numZplanes
            reg_Stacks{zStack} = reg___Stacks{zStack}(:,:,1:cutOffFrame);
%             CumData2{zStack}(:,:,cutOffFrame+1:end) = [];
%             CumDffDataFirst20s2{zStack}(:,:,cutOffFrame+1:end) = [];
%             dffDataFirst20s2{zStack}(:,:,cutOffFrame+1:end) = [];
        end        
    end 
    
    ResampedVel_wheel__data = ResampedVel_wheel_data; clear ResampedVel_wheel_data; 
    ResampedVel_wheel_data = ResampedVel_wheel__data(1:cutOffFrame);
end 

%% separate stacks by zPlane and trial type 
disp('Organizing Z-Stacks by Trial Type')
%go to the right directory for functions 
[imAn1funcDir] = getUserInput(userInput,'imAnalysis1_functions Directory');
cd(imAn1funcDir);

%find the diffent trial types 
[stimTimes] = getUserInput(userInput,"Stim Time Lengths (sec)"); 
[stimTypeNum] = getUserInput(userInput,"How many different kinds of stimuli were used?");
[uniqueTrialData,uniqueTrialDataOcurr,indices,state_start_f,~] = separateTrialTypes(TrialTypes,state_start_f,state_end_f,stimTimes,numZplanes,FPS,stimTypeNum);

if volIm == 1
    %separate the Z-stacks 
    sortedStacks2 = cell(1,length(reg_Stacks));
    for Zstack = 1:length(reg_Stacks)
          [sorted_Stacks,indices] = eventTriggeredAverages_STACKS2(reg_Stacks{Zstack},state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
          sortedStacks2{Zstack} = sorted_Stacks;           
    end 
elseif volIm == 0
    %separate the Z-stacks     
      [sorted_Stacks,indices] = eventTriggeredAverages_STACKS2(reg_Stacks{1},state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
      sortedStacks2{1} = sorted_Stacks;           
end 

[sortedStacks2,~,~] = removeEmptyCells(sortedStacks2,indices);

%below removes indices in cells where sortedStacks is blank 
for trialType = 1:size(sortedStacks2{1},2)
    if isempty(sortedStacks2{1}{trialType}) == 1 
        indices{trialType} = [];         
    end 
    if size(sortedStacks2{1}{trialType},2) < size(indices{trialType},1)
        indices{trialType}(size(sortedStacks2{1}{trialType},2)+1:end,:) = [];
    end 
end 

if cumStacksQ == 1  
    [dffStacks2,CumDffStacks2,CumStacks2] = makeCumPixStacksPerTrial(sortedStacks2,FPS,numZplanes,sec_before_stim_start);
end 


%% vessel segmentation 
if VsegQ == 1 
   [sortedData,userInput] = segmentVessels2(reg_Stacks,userInput,state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,numZplanes,ROIboundData);
end 

%% measure changes in calcium dynamics and BBB permeability 
if pixIntQ == 1
    if CaQ == 1
        %find max number of cells/terminals 
        maxCells = length(ROIinds);
        %determine change in pixel intensity sorted by cell identity
        %across Z 
        meanPixIntArray = cell(1,ROIinds(maxCells));
        for ccell = 1:maxCells %cell starts at 2 because that's where our cell identity labels begins (see identifyROIsAcrossZ function)
            %find the number of z planes a cell/terminal appears in 
            count = 1;
            %this figures out what planes in Z each cell occurs in (cellZ)
            for Z = 1:length(CaROImasks)                
                if ismember(ROIinds(ccell),CaROImasks{Z}) == 1 
                    cellInd = max(unique(ROIorders{Z}(CaROImasks{Z} == ROIinds(ccell))));
                    for frame = 1:length(reg_Stacks{Z})
                        stats = regionprops(ROIorders{Z},reg_Stacks{Z}(:,:,frame),'MeanIntensity');
                        meanPixIntArray{ROIinds(ccell)}(Z,frame) = stats(cellInd).MeanIntensity;
                    end 
                end 
            end 
            %turn all rows of zeros into NaNs
            allZeroRows = find(all(meanPixIntArray{ROIinds(ccell)} == 0,2));
            for row = 1:length(allZeroRows)
                meanPixIntArray{ROIinds(ccell)}(allZeroRows(row),:) = NaN; 
            end 
        end 

         dataMeds = cell(1,ROIinds(maxCells));
         DFOF = cell(1,ROIinds(maxCells));
         dataSlidingBLs = cell(1,ROIinds(maxCells));
         DFOFsubSBLs = cell(1,ROIinds(maxCells));
         zData = cell(1,ROIinds(maxCells));
         for ccell = 1:maxCells     
                for z = 1:size(meanPixIntArray{ROIinds(ccell)},1)     
                    %get median value per trace
                    dataMed = median(meanPixIntArray{ROIinds(ccell)}(z,:));     
                    dataMeds{ROIinds(ccell)}(z,:) = dataMed;
                    %compute DF/F using means  
                    DFOF{ROIinds(ccell)}(z,:) = (meanPixIntArray{ROIinds(ccell)}(z,:)-dataMeds{ROIinds(ccell)}(z,:))./dataMeds{ROIinds(ccell)}(z,:);                         
                    %get sliding baseline 
                    [dataSlidingBL]=slidingBaseline(DFOF{ROIinds(ccell)}(z,:),floor((FPS/numZplanes)*10),0.5); %0.5 quantile thresh = the median value                 
                    dataSlidingBLs{ROIinds(ccell)}(z,:) = dataSlidingBL;                       
                    %subtract sliding baseline from DF/F
                    DFOFsubSBLs{ROIinds(ccell)}(z,:) = DFOF{ROIinds(ccell)}(z,:)-dataSlidingBLs{ROIinds(ccell)}(z,:);
                    %z-score data 
                    zData{ROIinds(ccell)}(z,:) = zscore(DFOFsubSBLs{ROIinds(ccell)}(z,:));
                end
         end        

        %sort calcium data by trial type 
        for ccell = 1:maxCells
            for Z = 1:size(zData{ROIinds(ccell)},1)                
                [sortedStatArray,indices] = eventTriggeredAverages(zData{ROIinds(ccell)}(Z,:),state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);            
                sortedData{ROIinds(ccell)}(Z,:) = sortedStatArray;
            end                  
        end    
    end    
end 


%% wheel data work goes here 
%get median wheel value 
WdataMed = median(ResampedVel_wheel_data);     
%compute Dv/v using means  
DVOV = (ResampedVel_wheel_data-WdataMed)./WdataMed;      
%get sliding baseline 
[WdataSlidingBL]=slidingBaseline(DVOV,floor((FPS/numZplanes)*10),0.5); %0.5 quantile thresh = the median value                                   
%subtract sliding baseline from Dv/v
DVOVsubSBLs = DVOV-WdataSlidingBL;
%z-score wheel data 
zWData = zscore(DVOVsubSBLs);
%sort wheel data                    
[sortedWheelData,~] = eventTriggeredAverages(zWData,state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
 
 %% resort data by trial type 
 
 %resort the indices 
 indS = cell(1,size(sortedStacks2{1},2));
 indI = cell(1,size(sortedStacks2{1},2));
 for trialType = 1:size(sortedStacks2{1},2)   
    [S, I] = sort(indices{trialType});
    indS{trialType} = S;
    indI{trialType} = I;
 end 
 
%figure out max column of data types available 
tTypeInds = zeros(1,size(sortedStacks2{1},2));
for trialType = 1:size(sortedStacks2{1},2)
    if ismember(uniqueTrialData(trialType,:),uniqueTrialDataTemplate,'rows') == 1 
        [~, idx] = ismember(uniqueTrialData(trialType,:),uniqueTrialDataTemplate,'rows');
        tTypeInds(trialType) = idx; 
    end 
end 
maxTtypeInd = max(tTypeInds);

%sort data into correct based on trial type
if pixIntQ == 1 
    sortedData2 = cell(1,ROIinds(maxCells));
    for ccell = 1:maxCells    
        for z = 1:size(sortedData{ROIinds(ccell)},1)  
            for trialType = 1:maxTtypeInd 
                if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
                    [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
                    [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
                    
                    indI2{idxFin} = indI{idxStart};
                    indices2{idxFin} = indices{idxStart};                   
                    sortedData2{ROIinds(ccell)}{z,idxFin} = sortedData{ROIinds(ccell)}{z,idxStart};                    
                end 
            end
        end 
    end 
    indices2 = indices2';
end 

if VsegQ == 1
    for z = 1:length(sortedData)
        for ROI = 1:size(sortedData{1},2)
            for trialType = 1:maxTtypeInd        
                if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
                    [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
                    [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
                    
                    indI2{idxFin} = indI{idxStart};
                    indices2{idxFin} = indices{idxStart};                   
                    sortedData2{z}{ROI}{idxFin} = sortedData{z}{ROI}{idxStart};                    
                end 
            end          
        end 
    end 
    indices2 = indices2';
end 

if cumStacksQ == 1
    sortedStacks3 = cell(1,length(sortedStacks));
%     CumStacks3 = cell(1,length(sortedStacks));
%     CumDffStacks3 = cell(1,length(sortedStacks));
%     dffStacks3 = cell(1,length(sortedStacks));
    for z = 1:length(sortedStacks)
        for trialType = 1:maxTtypeInd     
            if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
                [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
                [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');

                indI2{idxFin} = indI{idxStart};
                indices2{idxFin} = indices{idxStart};                   
                sortedStacks3{z}{idxFin} = sortedStacks2{z}{idxStart};  
%                 CumStacks3{z}{idxFin} = CumStacks2{z}{idxStart};
%                 CumDffStacks3{z}{idxFin} = CumDffStacks2{z}{idxStart};
%                 dffStacks3{z}{idxFin} = dffStacks2{z}{idxStart};
            end 
        end          
    end 
    indices2 = indices2';
end 

%sort wheel data into correct spot based on trial type 
for trialType = 1:maxTtypeInd
    if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
        [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
        [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
                  
        sortedWheelData2{idxFin} = sortedWheelData{idxStart};                    
    end     
end

%% reorganize data by trial order

if pixIntQ == 1     
    dataToPlot2 = cell(1,ROIinds(maxCells));
     for ccell = 1:maxCells
        for z = 1:size(sortedData2{ROIinds(ccell)},1)
            for trialType = 1:size(indices2,1)
                if isempty(sortedData2{ROIinds(ccell)}{z,trialType}) == 0
                    dataToPlot2{ROIinds(ccell)}{z,trialType} = sortedData2{ROIinds(ccell)}{z,trialType}(indI2{trialType});       
                end 
            end 
        end 
     end 
end 

if VsegQ == 1
    %reorganize data by trial order 
    for z = 1:length(sortedData)
        for ROI = 1:size(sortedData{1},2)
            for trialType = 1:size(indices2,1)
                if isempty(sortedData2{z}{ROI}{trialType}) == 0
                    dataToPlot2{z}{ROI}{trialType} = sortedData2{z}{ROI}{trialType}(indI2{trialType});       
                end 
            end 
        end 
     end 
end 

if cumStacksQ == 1 
    clear sortedStacks2 CumStacks2 CumDffStacks2 dffStacks2;
    sortedStacks2 = cell(1,length(sortedStacks));
%     CumStacks2 = cell(1,length(sortedStacks));
%     CumDffStacks2 = cell(1,length(sortedStacks));
%     dffStacks2 = cell(1,length(sortedStacks));
    for z = 1:length(sortedStacks)        
        for trialType = 1:maxTtypeInd     
            if isempty(indices2{trialType}) == 0
                sortedStacks2{z}{trialType} = sortedStacks3{z}{trialType}(indI2{trialType});    
%                 CumStacks2{z}{trialType} = CumStacks3{z}{trialType}(indI2{trialType}); 
%                 CumDffStacks2{z}{trialType} = CumDffStacks3{z}{trialType}(indI2{trialType}); 
%                 dffStacks2{z}{trialType} = dffStacks3{z}{trialType}(indI2{trialType}); 
            end   
        end  
    end 
end 

wheelDataToPlot2 = cell(1,maxTtypeInd);
for trialType = 1:size(indices2,1)
    if isempty(indices2{trialType}) == 0 
        wheelDataToPlot2{trialType} = sortedWheelData2{trialType}(indI2{trialType});     
    end
end 

%% concatenate data 

if pixIntQ == 1 
    [dataToPlot3] = catData(dataToPlot,dataToPlot2,maxCells,ROIinds);
end 

if VsegQ == 1
    if size(dataToPlot2{1}{1},2) < size(dataToPlot{1}{1},2)
        ind1 = size(dataToPlot2{1}{1},2) + 1; ind2 = size(dataToPlot{1}{1},2);
        for z = 1:length(dataToPlot)
            for ROI = 1:size(dataToPlot{1},2)
                dataToPlot2{z}{ROI}{ind1:ind2} = [];
            end 
        end 
    end 
    
    [dataToPlot3] = catData2(dataToPlot,dataToPlot2);
end 

if cumStacksQ == 1
    clear CumStacks3 CumDffStacks3 dffStacks3 sortedStacks3
%     [CumStacks3] = catData3(CumStacks,CumStacks2);
%     [CumDffStacks3] = catData3(CumDffStacks,CumDffStacks2);
%     [dffStacks3] = catData3(dffStacks,dffStacks2);
    [sortedStacks3] = catData3(sortedStacks,sortedStacks2);
%     CumData3 = vertcat(CumData,CumData2);
%     CumDffDataFirst20s3 = vertcat(CumDffDataFirst20s,CumDffDataFirst20s2);
%     dffDataFirst20s3 = vertcat(dffDataFirst20s,dffDataFirst20s2);
end 

[wheelDataToPlot3] = catWheelData(wheelDataToPlot,wheelDataToPlot2);
 

%% prep data for plotting - get rid of what you don't need and average 
if pixIntQ == 1 || VsegQ == 1
    dataToPlot = dataToPlot3;
    clear dataToPlot3 
end 
        
wheelDataToPlot = wheelDataToPlot3; 
clear wheelDataToPlot3

if pixIntQ == 1    
    AVsortedData = cell(1,ROIinds(maxCells));
    for ccell = 1:maxCells  
        for z = 1:size(dataToPlot{ROIinds(ccell)},1)
            for trialType = 1:size(dataToPlot{ROIinds(ccell)},2) 
                if isempty(dataToPlot{ROIinds(ccell)}{z,trialType}) == 0 
                    reshapedArray = cat(3,dataToPlot{ROIinds(ccell)}{z,trialType}{:});
                    AV = nanmean(reshapedArray,3);
                    AVsortedData{ROIinds(ccell)}{z,trialType}(1,:) = AV; 
                end 
            end 
        end 
    end        
end 

if VsegQ == 1   
    for z = 1:length(sortedData)
        for ROI = 1:size(sortedData{1},2)
            for trialType = 1:size(dataToPlot{1}{1},2) 
                 if isempty(dataToPlot{z}{ROI}{trialType}) == 0 
                    reshapedArray = cat(3,dataToPlot{z}{ROI}{trialType}{:});
                    AV = nanmean(reshapedArray,3);
                    AVsortedData{z}{ROI}{trialType}(1,:) = AV; 
                 end 
 
            end 
        end 
    end      
    maxCells = size(sortedData{1},2);
end 

% average dff, cum dff, and cum stacks across all trials ALSO average pix intensity/vessel width data if applicable
if cumStacksQ == 1 
    clear sorted_Stacks
%     clear AVcumDffStacks AVcumStacks AVdffStacks AVStacks
%     CumDff_Stacks = cell(1,numZplanes);
%     Cum_Stacks = cell(1,numZplanes);
%     dff_Stacks = cell(1,numZplanes);
%     AVcumDffStacks = cell(1,numZplanes);
%     AVcumStacks = cell(1,numZplanes);
%     AVdffStacks = cell(1,numZplanes);
    AVStacks = cell(1,numZplanes);
    for Z = 1:numZplanes
        for trialType = 1:size(sortedStacks3{Z},2)
            if isempty(sortedStacks3{Z}{trialType}) == 0 
                for trial = 1:length(sortedStacks3{Z}{trialType})
%                     CumDff_Stacks{Z}{trialType}(:,:,:,trial) = CumDffStacks3{Z}{trialType}{trial};
%                     Cum_Stacks{Z}{trialType}(:,:,:,trial) = CumStacks3{Z}{trialType}{trial};
%                     dff_Stacks{Z}{trialType}(:,:,:,trial) = dffStacks3{Z}{trialType}{trial};
                    sorted_Stacks{Z}{trialType}(:,:,:,trial) = sortedStacks3{Z}{trialType}{trial};
                end 
%                 AVcumDffStacks{Z}{trialType} = mean(CumDff_Stacks{Z}{trialType},4);
%                 AVcumStacks{Z}{trialType} = mean(Cum_Stacks{Z}{trialType},4);
%                 AVdffStacks{Z}{trialType} = mean(dff_Stacks{Z}{trialType},4);
                AVStacks{Z}{trialType} = mean(sorted_Stacks{Z}{trialType},4);
            end 
        end 
    end 
end 

AVwheelData = cell(1,length(wheelDataToPlot));
for trialType = 1:length(wheelDataToPlot)
    if isempty(wheelDataToPlot{trialType}) == 0 
        reshapedArray = cat(3,wheelDataToPlot{trialType}{:});
        AV = nanmean(reshapedArray,3);
        AVwheelData{trialType}(1,:) = AV;  
    end 
end 

 %% clear unecessary values 

if cumStacksQ == 1
    %CumDffStacks = CumDffStacks3; CumStacks = CumStacks3; dffStacks = dffStacks3; sortedStacks = sortedStacks3; CumData = CumData3; CumDffDataFirst20s = CumDffDataFirst20s3; dffDataFirst20s = dffDataFirst20s3;
    sortedStacks = sortedStacks3; 
    clearvars -except dataToPlot AVsortedData wheelDataToPlot AVwheelData userInput FPS dataMin dataMax velMin velMax HDFchart numZplanes BG_ROIboundData CaROImasks uniqueTrialDataTemplate maxCells ROIorders ROIinds ROIboundData sec_before_stim_start sortedStacks CumStacks CumDffStacks dffStacks CumData CumDffDataFirst20s dffDataFirst20s AVcumDffStacks AVcumStacks AVdffStacks AVStacks
elseif  VsegQ == 1 || pixIntQ == 1
    clearvars -except dataToPlot AVsortedData wheelDataToPlot AVwheelData userInput FPS dataMin dataMax velMin velMax HDFchart numZplanes BG_ROIboundData CaROImasks uniqueTrialDataTemplate maxCells ROIorders ROIinds ROIboundData sec_before_stim_start
end 

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%THE BELOW CODE IS FOR PLOTTING - HAS NOT BEEN STREAMLINED. LOOK AT IT AND
%USE WHAT YOU NEED. 

% dataMin = input("data Y axis MIN: ");
% dataMax = input("data Y axis MAX: ");
% velMin = input("running velocity Y axis MIN: ");
% velMax = input("running velocity Y axis MAX: ");

%V = input('What vessel (diameter) do you want to see? '); 

%plotDataAndRunVelocity(dataToPlot,AVsortedData,wheelDataToPlot,AVwheelData,FPS,numZplanes,sec_before_stim_start,dataMin,dataMax,velMin,velMax,maxCells,ROIinds)
%plotAVDataAndRunVelocity(VdataToPlot,VAVsortedData,dataToPlot,AVsortedData,wheelDataToPlot,AVwheelData,FPS,numZplanes,sec_before_stim_start,dataMin,dataMax,maxCells,ROIinds,V)
%plotAVtTypeDataAndRunVelocity(VdataToPlot,VAVsortedData,AVtType1,AVtType2,AVtType3,AVtType4,AVAVtType1,AVAVtType2,AVAVtType3,AVAVtType4,wheelDataToPlot,AVwheelData,FPS,numZplanes,sec_before_stim_start,dataMin,dataMax,velMin,velMax,maxCells,ROIinds)
%plotAllAVDataAndRunVelocity(VdataToPlot,VAVsortedData,allAVarray,allAV,wheelDataToPlot,AVwheelData,FPS,numZplanes,sec_before_stim_start,dataMin,dataMax,velMin,velMax,maxCells,ROIinds)

%% in case the V data array is missing some trial type data and it needs to be the same size as the calcium data array for plotting 
%{
SEMVdata = cell(1,length(VAVsortedData));
for z = 1:length(VAVsortedData)
        for ROI = 1:size(VAVsortedData{z},2)
            for trialType = 1:size(AVsortedData{2},2)
                if size(VAVsortedData{z}{ROI},2) < trialType
                    for frame = 1:635
                        %VAVsortedData{z}{ROI}{trialType}(1,frame) = NaN; 
                        SEMVdata{z}{ROI}{trialType}(1,frame) = NaN;
                    end                    
                end 
            end 
        end 
end 

% 
%% average by trial types  
[AVsortedData2] = makeCellArrayCellsSameSize(AVsortedData,maxCells,ROIinds);

%average all trial types across all ROIs (calcium data) 
terminals = [15,13,11,10,2,3];

for ccell = 1:length(terminals)  
    for z = 1%:3        
        tType1(z,:,ccell) = AVsortedData2{terminals(ccell)}{z,1}; 
%         tType2(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,2}; 
%         tType3(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,3}; 
        %tType4(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,4}; 
    end 
 end 

% for ccell = 1:maxCells  
%     for z = 1%:3        
%         tType1(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,1}; 
% %         tType2(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,2}; 
% %         tType3(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,3}; 
%         %tType4(z,:,ccell) = AVsortedData2{ROIinds(ccell)}{z,4}; 
%     end 
%  end 
 
AVtType1 = nanmean(tType1,3); 
AVtType2 = nanmean(tType2,3);  
AVtType3 = nanmean(tType3,3); 
%AVtType4 = nanmean(tType4,3); 

AVAVtType1 = nanmean(AVtType1,1); 
AVAVtType2 = nanmean(AVtType2,1);  
AVAVtType3 = nanmean(AVtType3,1); 
%AVAVtType4 = nanmean(AVtType4,1); 



allAVarray(1,:) = AVAVtType1; %allAVarray(2,:) = AVAVtType2(1,1:size(AVsortedData{2}{1,1},2)); allAVarray(3,:) = AVAVtType3; %allAVarray(4,:) = AVAVtType4(1,1:size(AVsortedData{2}{1,1},2));
allAV = mean(allAVarray,1);

%per terminal - average all the trials 

cellAVtType1 = nanmean(tType1,1);
cellAVtType2 = nanmean(tType2,1);
cellAVtType3 = nanmean(tType3,1);
%cellAVtType4 = nanmean(tType4,1);

cellAllAv(1,:,:) = cellAVtType1; %cellAllAv(2,:,:) = cellAVtType2(:,1:size(AVsortedData{2}{1,1},2),:); cellAllAv(3,:,:) = cellAVtType3;%cellAllAv(4,:,:) = cellAVtType4(:,1:size(AVsortedData{2}{1,1},2),:);
cellAllAvRed(1,:,:) = cellAVtType3;%cellAllAvRed(2,:,:) = cellAVtType4(:,1:size(AVsortedData{2}{1,1},2),:);
cellAllAvBlue(1,:,:) = cellAVtType1;cellAllAvBlue(2,:,:) = cellAVtType2(:,1:size(AVsortedData{2}{1,1},2),:);

cellAllAvAV = nanmean(cellAllAv,1);
cellAllAvAVred = nanmean(cellAllAvRed,1);
cellAllAvAVblue = nanmean(cellAllAvBlue,1);

ALLRED = nanmean(cellAllAvAVred,3);
ALLBLUE = nanmean(cellAllAvAVblue,3);

%separate vessel width data by trial type - red vs blue 
%average all trial types across all ROIs (calcium data) 
V = 1;
redV(1,:) = VAVsortedData{1}{V}{3}; redV(2,:) = VAVsortedData{2}{V}{3}; redV(3,:) = VAVsortedData{3}{V}{3};
%redV(4,:) = VAVsortedData{1}{V}{4}(1,1:size(AVsortedData{2}{1,1},2)); redV(5,:) = VAVsortedData{2}{V}{4}(1,1:size(AVsortedData{2}{1,1},2)); redV(6,:) = VAVsortedData{3}{V}{4}(1,1:size(AVsortedData{2}{1,1},2));

V = 2;
blueV(1,:) = VAVsortedData{1}{V}{1}; blueV(2,:) = VAVsortedData{2}{V}{1}; blueV(3,:) = VAVsortedData{3}{V}{1};
blueV(4,:) = VAVsortedData{1}{V}{2}(1,1:size(AVsortedData{2}{1,1},2)); blueV(5,:) = VAVsortedData{2}{V}{2}(1,1:size(AVsortedData{2}{1,1},2)); blueV(6,:) = VAVsortedData{3}{V}{2}(1,1:size(AVsortedData{2}{1,1},2));


redVAV = nanmean(redV,1);
blueVAV = nanmean(blueV,1);

allVAVarray(1,:) = redVAV ; allVAVarray(2,:) = blueVAV;
allVAV = nanmean(allVAVarray,1);

for trial = 1:length(VdataToPlot{1}{1}{1})
    VDataArray(:,:,trial) =VdataToPlot{1}{V}{1}{trial};
end 

%average wheel data: RED, BLUE, and REDBLUE 
redW(1,:) = AVwheelData{3}; %redW(2,:) = AVwheelData{4}(1:size(AVwheelData{1},2));
blueW(1,:) = AVwheelData{1}; blueW(2,:) = AVwheelData{2}(1:size(AVwheelData{1},2));

AVredW = nanmean(redW,1);
AVblueW = nanmean(blueW,1);

AllW(1,:) = AVredW; AllW(2,:) = AVblueW;
AVallW = nanmean(AllW,1);

for trial = 1:length(wheelDataToPlot{1})
    wheelDataArray(:,:,trial) = wheelDataToPlot{1}{trial};
end 

%}
