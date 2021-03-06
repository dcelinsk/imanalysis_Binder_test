% make sure we know where the functions are 
imAnalysisDir = uigetdir('*.*','WHERE IS THE imAnalysis FOLDER? ');  
imAn1str = '\imAnalysis1_functions';
imAn2str = '\imAnalysis2_functions';
imAnDirs = [imAnalysisDir,imAn1str;imAnalysisDir,imAn2str];

%% register images
[regStacks,userInput,UIr,state_start_f,state_end_f,vel_wheel_data,TrialTypes,HDFchart] = imRegistration(imAnDirs);

%% set what data you want to plot 
dataParseType = input('What data do you need? Peristimulus epoch = 0. Stimulus epoch = 1. '); userInput(UIr,1) = ("What data do you need? Peristimulus epoch = 0. Stimulus epoch = 1."); userInput(UIr,2) = (dataParseType); UIr = UIr+1;    
if dataParseType == 0 
    sec_before_stim_start = input("How many seconds before the stimulus starts do you want to plot? "); userInput(UIr,1) = ("How many seconds before the stimulus starts do you want to plot?"); userInput(UIr,2) = (sec_before_stim_start); UIr = UIr+1;
    sec_after_stim_end = input("How many seconds after stimulus end do you want to plot? "); userInput(UIr,1) = ("How many seconds after stimulus end do you want to plot?"); userInput(UIr,2) = (sec_after_stim_end); UIr = UIr+1;
end 

%% set up analysis pipeline 
VsegQ = input('Do you need to measure vessel width? Yes = 1. No = 0. '); userInput(UIr,1) = ("Do you need to measure vessel width? Yes = 1. No = 0."); userInput(UIr,2) = (VsegQ); UIr = UIr+1;    
pixIntQ = input('Do you need to measure changes in pixel intensity? Yes = 1. No = 0. '); userInput(UIr,1) = ("Do you need to measure changes in pixel intensity? Yes = 1. No = 0."); userInput(UIr,2) = (pixIntQ); UIr = UIr+1; 
if pixIntQ == 1
    CaQ = input('Do you need to measure changes in calcium dynamics? Yes = 1. No = 0. '); userInput(UIr,1) = ("Do you need to measure changes in calcium dynamics? Yes = 1. No = 0."); userInput(UIr,2) = (CaQ); UIr = UIr+1; 
    BBBQ = input('Do you need to measure changes in BBB permeability? Yes = 1. No = 0. '); userInput(UIr,1) = ("Do you need to measure changes in BBB permeability? Yes = 1. No = 0."); userInput(UIr,2) = (BBBQ); UIr = UIr+1; 
    %BBBQ is there in case I want to make ROIs for measuring change in
    %pixel intensity of the cumStacks - this code has not been added in yet
end 
cumStacksQ = input('Do you want to generate cumulative pixel intensity stacks? Yes = 1. No = 0. '); userInput(UIr,1) = ("Do you want to generate cumulative pixel intensity stacks? Yes = 1. No = 0."); userInput(UIr,2) = (cumStacksQ); UIr = UIr+1; 

%% select registration method that's most appropriate for making the dff and cum pix int stacks 
[volIm] = getUserInput(userInput,'Is this volume imaging data? Yes = 1. Not = 0.');
if cumStacksQ == 1 || pixIntQ == 1 
    if volIm == 0
        regTypeDim = 0; userInput(UIr,1) = ("What registration dimension is best for pixel intensity analysis? 2D = 0. 3D = 1."); userInput(UIr,2) = (regTypeDim); UIr = UIr+1;
    elseif volIm == 1 
        regTypeDim = input("What registration dimension is best for pixel intensity analysis? 2D = 0. 3D = 1. "); userInput(UIr,1) = ("What registration dimension is best for pixel intensity analysis? 2D = 0. 3D = 1."); userInput(UIr,2) = (regTypeDim); UIr = UIr+1;
    end 
    regTypeTemp = input("What registration template is best for pixel intensity analysis? red = 0. green = 1. "); userInput(UIr,1) = ("What registration template is best for pixel intensity analysis? red = 0. green = 1."); userInput(UIr,2) = (regTypeTemp); UIr = UIr+1;

    [reg__Stacks] = pickRegStack(regStacks,regTypeDim,regTypeTemp);
    [reg_Stacks,BG_ROIboundData] = backgroundSubtraction(reg__Stacks);
end
%% select registration method that's most appropriate for vessel segmentation 
if VsegQ == 1
    if volIm == 0
        regTypeDimVesSeg = 0; userInput(UIr,1) = ("What registration dimension is best for vessel segmentation? 2D = 0. 3D = 1."); userInput(UIr,2) = (regTypeDimVesSeg); UIr = UIr+1;
    elseif volIm == 1 
        regTypeDimVesSeg = input("What registration dimension is best for vessel segmentation? 2D = 0. 3D = 1. "); userInput(UIr,1) = ("What registration dimension is best for vessel segmentation? 2D = 0. 3D = 1."); userInput(UIr,2) = (regTypeDimVesSeg); UIr = UIr+1;
    end 
    regTypeTempVesSeg = input("What registration template is best for vessel segmentation? red = 0. green = 1. "); userInput(UIr,1) = ("What registration template is best for vessel segmentation? red = 0. green = 1."); userInput(UIr,2) = (regTypeTempVesSeg); UIr = UIr+1;

    [reg__StacksVesSeg] = pickRegStack(regStacks,regTypeDimVesSeg,regTypeTempVesSeg);
    [reg_Stacks,BG_ROIboundData] = backgroundSubtraction(reg__StacksVesSeg);
end 

%% make cumulative, diff-cumulative, and DF/F stacks to output for calcium and BBB perm analysis 
[numZplanes] = getUserInput(userInput,"How many Z planes are there?");
% [FPS] = getUserInput(userInput,"FPS"); 
% if cumStacksQ == 1  
%     [dffDataFirst20s,CumDffDataFirst20s,CumData] = makeCumPixWholeEXPstacks(FPS,reg_Stacks,numZplanes,sec_before_stim_start);
% end

%% make sure state start and end frames line up 
[state_start_f,state_end_f,TrialTypes] = makeSureStartEndTrialTypesLineUp(reg_Stacks,state_start_f,state_end_f,TrialTypes,numZplanes);

%% resample velocity data by trial type
if length(vel_wheel_data)*size(reg_Stacks{1},3) > 2^31
    ResampedVel_wheel_data1 = resample(vel_wheel_data,(round(length(vel_wheel_data)/8000)),length(vel_wheel_data));
    ResampedVel_wheel_data = resample(ResampedVel_wheel_data1,size(reg_Stacks{1},3),length(ResampedVel_wheel_data1));
elseif length(vel_wheel_data)*size(reg_Stacks{1},3) < 2^31
    ResampedVel_wheel_data = resample(vel_wheel_data,size(reg_Stacks{1},3),length(vel_wheel_data));
end 

%% get rid of frames/trials where registration gets wonky 
%EVENTUALLY MAKE THIS AUTOMATIC INSTEAD OF HAVING TO INPUT WHAT FRAME THE
%REGISTRATION GETS WONKY 
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
%             CumData{zStack}(:,:,cutOffFrame+1:end) = [];
%             CumDffDataFirst20s{zStack}(:,:,cutOffFrame+1:end) = [];
%             dffDataFirst20s{zStack}(:,:,cutOffFrame+1:end) = [];
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
[FPS] = getUserInput(userInput,"FPS");
[uniqueTrialData,uniqueTrialDataOcurr,indices,state_start_f,uniqueTrialDataTemplate] = separateTrialTypes(TrialTypes,state_start_f,state_end_f,stimTimes,numZplanes,FPS,stimTypeNum); 

if volIm == 1
    %separate the Z-stacks 
    sortedStacks = cell(1,length(reg_Stacks));
    for Zstack = 1:length(reg_Stacks)
          [sorted_Stacks,indices] = eventTriggeredAverages_STACKS2(reg_Stacks{Zstack},state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
          sortedStacks{Zstack} = sorted_Stacks;           
    end 
elseif volIm == 0
    %separate the Z-stacks      
      [sorted_Stacks,indices] = eventTriggeredAverages_STACKS2(reg_Stacks{1},state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
      sortedStacks{1} = sorted_Stacks;           
end 
        
[sortedStacks,~,~] = removeEmptyCells(sortedStacks,indices);

%below removes indices in cells where sortedStacks is blank 
for trialType = 1:size(sortedStacks{1},2)
    if isempty(sortedStacks{1}{trialType}) == 1 
        indices{trialType} = [];         
    end 
end 

if cumStacksQ == 1  
    [dffStacks,CumDffStacks,CumStacks] = makeCumPixStacksPerTrial(sortedStacks,FPS,numZplanes,sec_before_stim_start);
end

%% vessel segmentation
if VsegQ == 1
        [sortedData_width_Z,sortedData_width,sortedData_area,userInput,UIr,ROIboundData] = segmentVessels(reg_Stacks,volIm,UIr,userInput,state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,numZplanes,HDFchart);
end 

%% measure changes in calcium dynamics and BBB permeability 
% Dmitrijs
% if pixIntQ == 1
%     if CaQ == 1    
%         if volIm == 1
%             [CaROImasks,userInput,ROIorders] = identifyROIsAcrossZ(reg_Stacks,userInput,UIr,numZplanes);
%         elseif volIm == 0 
%             [CaROImasks,userInput,ROIorders] = identifyROIs(reg_Stacks,userInput,UIr);
%         end 
%         
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% %         %EDIT CALCIUM ROI MASKS BY HAND!!!!!
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% % %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% % 
% %           CaROImasks{1}(CaROImasks{1}==9)=0; CaROImasks{3}(CaROImasks{3}==11)=0; CaROImasks{1}(CaROImasks{1}==3)=0; CaROImasks{3}(CaROImasks{3}==8)=0;%CaROImasks{2}(CaROImasks{2}==10)=0;
% %           figure;imagesc(CaROImasks{1});grid on;figure;imagesc(CaROImasks{2});grid on;figure;imagesc(CaROImasks{3});grid on
% %         
%         masksDoneQ = input('Have the calcium ROI masks been hand edited? Yes = 1. No = 0.');
%         if masksDoneQ == 1 
%             %determine the indices left for the edited CaROImasks or else
%             %there will be indexing problems below through iteration 
%             ROIinds = unique([CaROImasks{:}]);
%             %remove zero
%             ROIinds(ROIinds==0) = [];
%             %find max number of cells/terminals 
%             maxCells = length(ROIinds);
%             %determine change in pixel intensity sorted by cell identity
%             %across Z 
%             meanPixIntArray = cell(1,ROIinds(maxCells));
%             for ccell = 1:maxCells %cell starts at 2 because that's where our cell identity labels begins (see identifyROIsAcrossZ function)
%                 %find the number of z planes a cell/terminal appears in 
%                 count = 1;
%                 %this figures out what planes in Z each cell occurs in (cellZ)
%                 for Z = 1:length(CaROImasks)                
%                     if ismember(ROIinds(ccell),CaROImasks{Z}) == 1 
%                         cellInd = max(unique(ROIorders{Z}(CaROImasks{Z} == ROIinds(ccell))));
%                         for frame = 1:length(reg_Stacks{Z})
%                             stats = regionprops(ROIorders{Z},reg_Stacks{Z}(:,:,frame),'MeanIntensity');
%                             meanPixIntArray{ROIinds(ccell)}(Z,frame) = stats(cellInd).MeanIntensity;
%                         end 
%                     end 
%                 end 
%                 %turn all rows of zeros into NaNs
%                 allZeroRows = find(all(meanPixIntArray{ROIinds(ccell)} == 0,2));
%                 for row = 1:length(allZeroRows)
%                     meanPixIntArray{ROIinds(ccell)}(allZeroRows(row),:) = NaN; 
%                 end 
%             end 
%             
%              dataMeds = cell(1,ROIinds(maxCells));
%              DFOF = cell(1,ROIinds(maxCells));
%              dataSlidingBLs = cell(1,ROIinds(maxCells));
%              DFOFsubSBLs = cell(1,ROIinds(maxCells));
%              zData = cell(1,ROIinds(maxCells));
%              for ccell = 1:maxCells     
%                     for z = 1:size(meanPixIntArray{ROIinds(ccell)},1)     
%                         %get median value per trace
%                         dataMed = nanmedian(meanPixIntArray{ROIinds(ccell)}(z,:));     
%                         dataMeds{ROIinds(ccell)}(z,:) = dataMed;
%                         %compute DF/F using means  
%                         DFOF{ROIinds(ccell)}(z,:) = (meanPixIntArray{ROIinds(ccell)}(z,:)-dataMeds{ROIinds(ccell)}(z,:))./dataMeds{ROIinds(ccell)}(z,:);                         
%                         %get sliding baseline 
%                         [dataSlidingBL]=slidingBaseline(DFOF{ROIinds(ccell)}(z,:),floor((FPS/numZplanes)*10),0.5); %0.5 quantile thresh = the median value                 
%                         dataSlidingBLs{ROIinds(ccell)}(z,:) = dataSlidingBL;                       
%                         %subtract sliding baseline from DF/F
%                         DFOFsubSBLs{ROIinds(ccell)}(z,:) = DFOF{ROIinds(ccell)}(z,:)-dataSlidingBLs{ROIinds(ccell)}(z,:);
%                         %z-score data 
%                         zData{ROIinds(ccell)}(z,:) = zscore(DFOFsubSBLs{ROIinds(ccell)}(z,:));
%                     end
%              end
% 
%             %sort calcium data by trial type 
%             for ccell = 1:maxCells
%                 for Z = 1:size(zData{ROIinds(ccell)},1)  
%                     [sortedStatArray,indices] = eventTriggeredAverages(zData{ROIinds(ccell)}(Z,:),state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);            
%                     sortedData{ROIinds(ccell)}(Z,:) = sortedStatArray;
%                 end                  
%             end            
%         end 
%     end 
% end

%% wheel data work goes here 
% %get median wheel value 
% WdataMed = median(ResampedVel_wheel_data);     
% %compute Dv/v using means  
% DVOV = (ResampedVel_wheel_data-WdataMed)./WdataMed;      
% %get sliding baseline 
% [WdataSlidingBL]=slidingBaseline(DVOV,floor((FPS/numZplanes)*10),0.5); %0.5 quantile thresh = the median value                                   
% %subtract sliding baseline from Dv/v
% DVOVsubSBLs = DVOV-WdataSlidingBL;
% %z-score wheel data 
% zWData = zscore(DVOVsubSBLs);
% %sort wheel data                    
% [sortedWheelData,~] = eventTriggeredAverages(zWData,state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);

%% average dff, cum dff,cum stacks across all trials, ALSO average pix intensity/vessel width data if applicable
% Dmitrijs
% if cumStacksQ == 1 
% %     CumDff_Stacks = cell(1,numZplanes);
% %     Cum_Stacks = cell(1,numZplanes);
% %     dff_Stacks = cell(1,numZplanes);
% %     AVcumDffStacks = cell(1,numZplanes);
% %     AVcumStacks = cell(1,numZplanes);
% %     AVdffStacks = cell(1,numZplanes);
%     AVStacks = cell(1,numZplanes);
%     for Z = 1:numZplanes
%         for trialType = 1:size(sortedStacks{1},2)
%             if isempty(sortedStacks{Z}{trialType}) == 0
%                 for trial = 1:length(CumDffStacks{Z}{trialType})
% % %                     CumDff_Stacks{Z}{trialType}(:,:,:,trial) = CumDffStacks{Z}{trialType}{trial};
% % %                     Cum_Stacks{Z}{trialType}(:,:,:,trial) = CumStacks{Z}{trialType}{trial};
% % %                     dff_Stacks{Z}{trialType}(:,:,:,trial) = dffStacks{Z}{trialType}{trial};
%                     sorted_Stacks{Z}{trialType}(:,:,:,trial) = sortedStacks{Z}{trialType}{trial};
%                 end 
% %                 AVcumDffStacks{Z}{trialType} = mean(CumDff_Stacks{Z}{trialType},4);
% %                 AVcumStacks{Z}{trialType} = mean(Cum_Stacks{Z}{trialType},4);
% %                 AVdffStacks{Z}{trialType} = mean(dff_Stacks{Z}{trialType},4);
%                 AVStacks{Z}{trialType} = mean(sorted_Stacks{Z}{trialType},4);
%             end 
%         
%         end 
%     end 
% end 
% 
% if pixIntQ == 1        
%     sortedStats_Array = cell(1,ROIinds(maxCells));
%     AVsortedData = cell(1,ROIinds(maxCells));
%     for ccell = 1:maxCells       
%         for z = 1:size(sortedData{ROIinds(ccell)},1)
%             for trialType = 1:size(sortedData{ROIinds(ccell)},2) 
%                 if isempty(sortedData{ROIinds(ccell)}{z,trialType}) == 0
%                     for trial = 1:length(sortedData{ROIinds(ccell)}{z,trialType})
%                          sortedStats_Array{ROIinds(ccell)}{z,trialType}(:,:,:,trial) = sortedData{ROIinds(ccell)}{z,trialType}{trial};
%                     end 
%                     AVsortedData{ROIinds(ccell)}{z,trialType}(1,:) = mean(sortedStats_Array{ROIinds(ccell)}{z,trialType},4);
%                 end 
%             end 
%         end 
%     end        
% end 

if VsegQ == 1
    for z = 1:length(sortedData_width_Z)
        for ROI = 1:size(sortedData_width_Z{1},2)
            for trialType = 1:size(sortedData_width_Z{1}{1},2)   
                if isempty(sortedData_width_Z{z}{ROI}{trialType}) == 0                  
                    for trial = 1:length(sortedData_width_Z{z}{ROI}{trialType})  
                        if isempty(sortedData_width_Z{z}{ROI}{trialType}{trial}) == 0  
                            sortedStats_Array{z}{ROI}{trialType}(:,:,:,trial) = sortedData_width_Z{z}{ROI}{trialType}{trial};
                        end 
                    end 
                    AVsortedData{z}{ROI}{trialType}(1,:) = mean(sortedStats_Array{z}{ROI}{trialType},4);
                end               
            end            
        end 
    end 
end 
% Dmitrijs
% sortedWheel_Data = cell(1,size(sortedStacks{1},2));
% AVwheelData = cell(1,size(sortedStacks{1},2));
% for trialType = 1:size(sortedStacks{1},2)   
%     if isempty(sortedWheelData{trialType}) == 0
%         for trial = 1:length(sortedWheelData{trialType})
%             if isempty(sortedWheelData{trialType}{trial}) == 0
%                 sortedWheel_Data{trialType}(:,:,:,trial) = sortedWheelData{trialType}{trial};
%             end 
%         end   
%         AVwheelData{trialType}(1,:) = mean(sortedWheel_Data{trialType},4); 
%     end 
% end 


%% reorganize data so that the trials are grouped together by type

%resort the indices 
 indS = cell(1,size(sortedStacks{1},2));
 indI = cell(1,size(sortedStacks{1},2));
 for trialType = 1:size(sortedStacks{1},2)   
    [S, I] = sort(indices{trialType});
    indS{trialType} = S;
    indI{trialType} = I;
 end 

%figure out max column of data types available 
tTypeInds = zeros(1,size(sortedStacks{1},2));
for trialType = 1:size(sortedStacks{1},2)
    if ismember(uniqueTrialData(trialType,:),uniqueTrialDataTemplate,'rows') == 1 
        [~, idx] = ismember(uniqueTrialData(trialType,:),uniqueTrialDataTemplate,'rows');
        tTypeInds(trialType) = idx; 
    end 
end 
maxTtypeInd = max(tTypeInds);

%sort data into correct spot based on trial type 
if pixIntQ == 1 
    sortedData2 = cell(1,ROIinds(maxCells));
    for ccell = 1:maxCells    
        for z = 1:size(sortedData_width_Z{ROIinds(ccell)},1)  
            for trialType = 1:maxTtypeInd 
                if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
                    [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
                    [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
                    
                    indI2{idxFin} = indI{idxStart};
                    indices2{idxFin} = indices{idxStart};                   
                    sortedData2{ROIinds(ccell)}{z,idxFin} = sortedData_width_Z{ROIinds(ccell)}{z,idxStart};                    
                end 
            end
        end 
    end 
    indices2 = indices2';
end 

if VsegQ == 1
    for z = 1:length(sortedData_width_Z)
        for ROI = 1:size(sortedData_width_Z{1},2)
            for trialType = 1:maxTtypeInd        
                if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
                    [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
                    [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
                    
                    indI2{idxFin} = indI{idxStart};
                    indices2{idxFin} = indices{idxStart};                   
                    sortedData2{z}{ROI}{idxFin} = sortedData_width_Z{z}{ROI}{idxStart};                    
                end 
            end          
        end 
    end 
    indices2 = indices2';
end 

if cumStacksQ == 1
    sortedStacks2 = cell(1,length(sortedStacks));
%     CumStacks2 = cell(1,length(sortedStacks));
%     CumDffStacks2 = cell(1,length(sortedStacks));
%     dffStacks2 = cell(1,length(sortedStacks));
%     AVStacks2 = cell(1,length(sortedStacks));
%     AVcumStacks2 = cell(1,length(sortedStacks));
%     AVcumDffStacks2 = cell(1,length(sortedStacks));
%     AVdffStacks2 = cell(1,length(sortedStacks));
    for z = 1:length(sortedStacks)
        for trialType = 1:maxTtypeInd     
            if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
                [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
                [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');

                indI2{idxFin} = indI{idxStart};
                indices2{idxFin} = indices{idxStart};                   
                sortedStacks2{z}{idxFin} = sortedStacks{z}{idxStart};  
%                 CumStacks2{z}{idxFin} = CumStacks{z}{idxStart};
%                 CumDffStacks2{z}{idxFin} = CumDffStacks{z}{idxStart};
%                 dffStacks2{z}{idxFin} = dffStacks{z}{idxStart};
%                 AVStacks2{z}{idxFin} = AVStacks{z}{idxStart};  
%                 AVcumStacks2{z}{idxFin} = AVcumStacks{z}{idxStart};
%                 AVcumDffStacks2{z}{idxFin} = AVcumDffStacks{z}{idxStart};
%                 AVdffStacks2{z}{idxFin} = AVdffStacks{z}{idxStart};
            end 
        end          
    end 
    indices2 = indices2';
end 

% Dmitrijs
% %sort wheel data into correct spot based on trial type 
% for trialType = 1:maxTtypeInd 
%     if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
%         [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
%         [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
%                   
%         sortedWheelData2{idxFin} = sortedWheelData{idxStart};                    
%     end     
% end

%% reorder trials from earliest to latest occurance in time %
% Dmitrijs
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% %@@@@@@@@@@@@@@@@@@    FIX THIS LATER    @@@@@@@@@@@@@@@@@@@@@@@@@@@@
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% % 
% dataToPlot = sortedData2;
% % wheelDataToPlot = sortedWheelData2;
% 
% % sortedStacks = sortedStacks2; %Dmitrijs
% % wheelDataToPlot = sortedWheelData2; %Dmitrijs
% 
% 
% % if pixIntQ == 1    
% %     dataToPlot = cell(1,ROIinds(maxCells));
% %     for ccell = 1:maxCells  
% %         for z = 1:size(sortedData2{ROIinds(ccell)},1)
% %             for trialType = 1:maxTtypeInd 
% %                 if isempty(sortedData2{ROIinds(ccell)}{z,trialType}) == 0
% %                     dataToPlot{ROIinds(ccell)}{z,trialType} = sortedData2{ROIinds(ccell)}{z,trialType}(indI2{trialType}(1:length(sortedData2{ROIinds(ccell)}{z,trialType})));     
% %                 end 
% %             end 
% %         end 
% %     end 
% % end 
% % 
% if VsegQ == 1  
%     for z = 1:length(sortedData)
%         for ROI = 1:size(sortedData{1},2)             
%             for trialType = 1:maxTtypeInd     
%                 if isempty(sortedData2{z}{ROI}{trialType}) == 0
%                     dataToPlot{z}{ROI}{trialType} = sortedData2{z}{ROI}{trialType}(indI2{trialType});     
%                 end   
%             end 
%         end 
%     end 
% end 
% % 
% % if cumStacksQ == 1 
% %     clear sortedStacks CumStacks CumDffStacks dffStacks AVStacks AVcumStacks AVcumDffStacks AVdffStacks;
% % %     AVStacks = AVStacks2; AVcumStacks = AVcumStacks2; AVcumDffStacks = AVcumDffStacks2; AVdffStacks = AVdffStacks2;
% %     sortedStacks = cell(1,length(sortedStacks2));
% %     CumStacks = cell(1,length(sortedStacks2));
% %     CumDffStacks = cell(1,length(sortedStacks2));
% %     dffStacks = cell(1,length(sortedStacks2));
% %     for z = 1:length(sortedStacks2)        
% %         for trialType = 1:maxTtypeInd     
% %             if isempty(sortedStacks2{z}{trialType}) == 0
% %                 sortedStacks{z}{trialType} = sortedStacks2{z}{trialType}(indI2{trialType});    
% % %                 CumStacks{z}{trialType} = CumStacks2{z}{trialType}(indI2{trialType}); 
% % %                 CumDffStacks{z}{trialType} = CumDffStacks2{z}{trialType}(indI2{trialType}); 
% % %                 dffStacks{z}{trialType} = dffStacks2{z}{trialType}(indI2{trialType}); 
% %             end   
% %         end  
% %     end 
% % end 
% %  
% % wheelDataToPlot = cell(1,maxTtypeInd);
% % for trialType = 1:maxTtypeInd 
% %     if isempty(sortedWheelData2{trialType}) == 0 
% %         wheelDataToPlot{trialType} = sortedWheelData2{trialType}(indI2{trialType});  
% %     end 
% % end 


%% Save ROI contour and traces
clc
variables_to_plot = {'sortedData_width','sortedData_width_Z','sortedData_area'};
y_labels = {'Vessel width (pix)', 'z-scored width', 'Vessel area (pix^2)'};
save_names = {'_mean_width_trace','_mean_widthZ_trace','_mean_area_trace'};
stim_names = {'blue','red'};
angle = getUserInput(userInput,"ROI Rotation Angles");
for m=1:length(variables_to_plot)
    dataToPlot = eval(string(variables_to_plot(m)));
    for i=1:length(ROIboundData)
        figure
        subplot(length(dataToPlot{1}{1})+1,1,1)
        imshow(imrotate(mean(reg_Stacks{1},3),angle(i)),[0,1800]);
        rectangle('Position',cell2mat(ROIboundData{i}),'EdgeColor','r','LineWidth',1)
        text(ROIboundData{i}{1}-16,ROIboundData{i}{2}+4,string(i),'FontSize',14,'Color','r')
        title({strcat('ROI', num2str(i))})

        for jj = 1:length(dataToPlot{1}{1})
            subplot(length(dataToPlot{1}{1})+1,1,jj+1)
            temp = reshape(cell2mat(dataToPlot{1}{1,i}{jj}),size(dataToPlot{1}{1,i}{jj}{2},2), []);
            t = (1:length(temp))/FPS; % in mins
            [t_bstim] = getUserInput(userInput,"How many seconds before the stimulus starts do you want to plot?");
            [t_stim] = getUserInput(userInput,"Stim Time Lengths (sec)");
            t = t-t_bstim;
            std_dev = std(temp,0,2);
            curve1 = mean(temp,2) + std_dev;
            curve2 = mean(temp,2) - std_dev;
            t2 = [t, fliplr(t)];
            inBetween = [curve1', fliplr(curve2')];
            plot(t, temp, 'color', [0.0285 0.6137 0.8135], 'LineWidth', 0.5);
            hold on;
            fill(t2, inBetween,[1 .9 .9],'linestyle','none','facealpha',0.7);
            plot(t, mean(temp,2), 'r', 'LineWidth', 2);
        
            xlabel('Time (s)')
            ylabel(y_labels(m))
            box off
            set(gca,'FontSize',12,'LineWidth',1.2)
            % Adding an injection shading
            x = [0 t_stim t_stim 0];
            y = [min(ylim) min(ylim) max(ylim) max(ylim)];
            p=patch(x,y,'b','EdgeColor','none');
            set(p,'FaceAlpha',0.5, 'FaceColor',[0 0.447 0.741]);
            title({strcat(string(stim_names(jj)),'Stim, n = ',...
                num2str(size(dataToPlot{1}{1,i}{jj},2)-1),' trials')})
            clear t std_dev curve1 curve2 t2 temp
        
            % Here we preserve the size of the image when we save it.
            % https://dgleich.github.io/hq-matlab-figs/
            % Minimize white space in the figure
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - 2*ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        end

        set(gcf,'InvertHardcopy','on');
        set(gcf,'PaperUnits', 'inches');
        width = 11.2; height = 8.4;
        papersize = get(gcf, 'PaperSize');
        left = (papersize(1)- width)/2;
        bottom = (papersize(2)- height)/2;
        myfiguresize = [left, bottom, width, height];
        set(gcf,'PaperPosition', myfiguresize);
        print(strcat('ROI',num2str(i),string(save_names(m))),'-dpng','-r600');
    end
    clear dataToPlot
end

%% Stitch trials together into continuous traces with gaps
clc
variables_to_plot = {'sortedData_width','sortedData_width_Z','sortedData_area'};
y_labels = {'Vessel width (pix)', 'z-scored width', 'Vessel area (pix^2)'};
save_names = {'_width_trace','_widthZ_trace','_area_trace'};
stim_names = {'blue','red'};
d = designfilt('lowpassfir', ...
    'PassbandFrequency',0.15,'StopbandFrequency',0.2, ...
    'PassbandRipple',1,'StopbandAttenuation',60, ...
    'DesignMethod','equiripple'); % Matlab example of ECG filter
for m=1:length(variables_to_plot)
    dataToPlot = eval(string(variables_to_plot(m)));
    for jj = 1:length(dataToPlot{1}{1})
        for i=1:length(ROIboundData)
            figure

            temp = reshape(cell2mat(dataToPlot{1}{1,i}{jj}),...
                size(dataToPlot{1}{1,i}{jj}{2},2), []);
            title({strcat(string(stim_names(jj)), 'Stim, n = ', num2str(size(temp,2)),' trials'),...
                strcat('ROI', num2str(i))})
            t = (1:length(temp))/FPS; % in sec
            [t_bstim] = getUserInput(userInput,"How many seconds before the stimulus starts do you want to plot?");
            dt = 4; % gaps between plots - actual gap is stim_period+t_bstim
            % added a gap just to make clear that this is not continuous data
            [t_stim] = getUserInput(userInput,"Stim Time Lengths (sec)");
            kk = 1;
            dt2 = 0; % for some reason every row in subplot lags by 10sec, but
    %         the time is not real anyways because of repeats and gaps
            for k = 1:3:length(temp(1,:))
                subplot(ceil(length(temp(1,:))/3),1,kk)
                for j = 0:2
                    if (k+j<=length(temp(1,:)))
        %                 plot(t+t_bstim*2*(k+j-1)+dt*j, filtfilt(d,temp(:,k+j)),...
        %                     'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);
                        plot(t+(t_bstim*2)*(k+j-1)+dt*j+dt2, temp(:,k+j),...
                            'color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1);
                        hold on;
                            % Adding an injection shading
                        x = [t_bstim+t_bstim*2*(k+j-1)+dt*j+dt2 t_bstim+t_stim+t_bstim*2*(k+j-1)+dt*j+dt2...
                            t_bstim+t_stim+t_bstim*2*(k+j-1)+dt*j+dt2 t_bstim+t_bstim*2*(k+j-1)+dt*j+dt2];
                        y = [min(ylim) min(ylim) max(ylim) max(ylim)];
                        p=patch(x,y,'b','EdgeColor','none');
                        set(p,'FaceAlpha',0.5, 'FaceColor',[0 0.447 0.741]);
                        ylim([min(min(temp)) max(max(temp))])
                    else
                        break
                    end
                end
                if kk==ceil(length(temp(1,:))/3/2)
                    ylabel(y_labels(m))
                end
                kk = kk+1;
                dt2=dt2+10;
            end

            suptitle({strcat(string(stim_names(jj)), 'Stim, n = ',...
                num2str(size(temp,2)),' trials'),strcat('ROI', num2str(i))}) % requires Bioinformatics toolbox
            xlabel('Time (s)')

            clear t std_dev curve1 curve2 t2 temp

                % Here we preserve the size of the image when we save it.
            % https://dgleich.github.io/hq-matlab-figs/
            set(gcf,'InvertHardcopy','on');
            set(gcf,'PaperUnits', 'inches');
            width = 11.2; height = 8.4;
            papersize = get(gcf, 'PaperSize');
            left = (papersize(1)- width)/2;
            bottom = (papersize(2)- height)/2;
            myfiguresize = [left, bottom, width, height];
            set(gcf,'PaperPosition', myfiguresize);
            print(strcat('ROI',num2str(i),string(save_names(m)),...
                '_', string(stim_names(jj)), 'Stim'),'-dpng','-r300');
        end
    end
    clear dataToPlot
end

%% clear unecessary values

if cumStacksQ == 1
    clearvars -except dataToPlot AVsortedData wheelDataToPlot AVwheelData userInput FPS dataMin dataMax velMin velMax HDFchart numZplanes BG_ROIboundData CaROImasks uniqueTrialDataTemplate maxCells ROIorders ROIinds ROIboundData sec_before_stim_start sortedStacks CumStacks CumDffStacks dffStacks CumData CumDffDataFirst20s dffDataFirst20s AVcumDffStacks AVcumStacks AVdffStacks AVStacks
elseif  VsegQ == 1 || pixIntQ == 1
    clearvars -except sortedData_width_Z sortedData_width sortedData_area AVsortedData wheelDataToPlot AVwheelData userInput FPS dataMin dataMax velMin velMax HDFchart numZplanes BG_ROIboundData CaROImasks uniqueTrialDataTemplate maxCells ROIorders ROIinds ROIboundData sec_before_stim_start
end 
%}

%% Save an image with all labeled ROIs
clc
figure
imshow(mean(reg_Stacks{1},3),[0,1800]);
for i=1:length(ROIboundData)
    rectangle('Position',cell2mat(ROIboundData{i}),'EdgeColor','r','LineWidth',2)
    text(ROIboundData{i}{1}-10,ROIboundData{i}{2}+4,string(i),'FontSize',14,'Color','r')
    hold on
end

% Here we preserve the size of the image when we save it.
% https://dgleich.github.io/hq-matlab-figs/
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
width = 11.2; height = 8.4;
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('ROI_contours','-dpng','-r600');
disp('finished')

%% Plot mean vessel width traces - shades
clc
dataToPlot = sortedData_width;

for i=1:length(ROIboundData)
    figure

    temp = reshape(cell2mat(dataToPlot{1}{1,i}{1}),...
        size(dataToPlot{1}{1,i}{1}{2},2), []);
    t = (1:length(temp))/FPS; % in mins
    [t_bstim] = getUserInput(userInput,"How many seconds before the stimulus starts do you want to plot?");
    [t_stim] = getUserInput(userInput,"Stim Time Lengths (sec)");
    t = t-t_bstim;
    std_dev = std(temp,0,2);
    curve1 = mean(temp,2) + std_dev;
    curve2 = mean(temp,2) - std_dev;
    t2 = [t, fliplr(t)];
    inBetween = [curve1', fliplr(curve2')];
    plot(t, temp, 'color', [0.0285 0.6137 0.8135], 'LineWidth', 0.01);
    hold on;
    fill(t2, inBetween,[1 .9 .9],'linestyle','none','facealpha',0.7);
    plot(t, mean(temp,2), 'r', 'LineWidth', 2);
    title({strcat('Stimulation at 470nm, n = ', num2str(size(temp,2)),' trials'),...
        strcat('ROI', num2str(i))})

    % xlim([-4 10])
    xlabel('Time (s)')
    ylabel('z-score')
    box off
    set(gca,'FontSize',12,'LineWidth',1.2)
    % Adding an injection shading
    x = [0 t_stim t_stim 0];
    y = [min(ylim) min(ylim) max(ylim) max(ylim)];
    p=patch(x,y,'b','EdgeColor','none');
    set(p,'FaceAlpha',0.5, 'FaceColor',[0 0.447 0.741]);

    clear t std_dev curve1 curve2 t2 temp
    
        % Here we preserve the size of the image when we save it.
    % https://dgleich.github.io/hq-matlab-figs/
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    width = 11.2; height = 8.4;
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    print(strcat('ROI',num2str(i),'_mean_trace'),'-dpng','-r300');
end


%% Plot the ROI box and amplitude across the trials right after the stim
clc
dataToPlot = sortedData_width;
angle = getUserInput(userInput,"ROI Rotation Angles");
dt = 1; % Time interval post-stimulus during which to average the amplitude
% dataToPlot = sortedData2;
[t_bstim] = getUserInput(userInput,"How many seconds before the stimulus starts do you want to plot?");
[t_stim] = getUserInput(userInput,"Stim Time Lengths (sec)");

for i=1:length(ROIboundData)
    figure
    subplot(2,1,1)
    imshow(imrotate(mean(reg_Stacks{1},3),angle(i)),[0,1800]);
    rectangle('Position',cell2mat(ROIboundData{i}),'EdgeColor','r','LineWidth',1)
    text(ROIboundData{i}{1}-16,ROIboundData{i}{2}+4,string(i),'FontSize',14,'Color','r')
    title({strcat('Stimulation at 470nm, ROI', num2str(i))})
    
    subplot(2,1,2)
    temp = reshape(cell2mat(dataToPlot{1}{1,i}{1}),...
        size(dataToPlot{1}{1,i}{1}{2},2), []);
    t = (1:length(temp))/FPS; % in sec
    amp_post = [];
    for k = 1:length(temp(1,:))
        % ind1...ind2 define the range from which we take the mean
        % amplitude
        [minValue1,ind1] = min(abs(t-(t_bstim+t_stim)));
        [minValue2,ind2] = min(abs(t-(t_bstim+t_stim+dt)));
%         t(ind2)
        amp_post(k) = mean(temp(ind1:ind2,k));
    end
    
    plot(1:length(temp(1,:)), amp_post,'--o',...
        'color', [0.9 0 0], 'LineWidth', 1, 'MarkerFaceColor',[0.9 0 0]);
    hold on
    plot(1:length(temp(1,:)), zeros(length(temp(1,:)),1),'--k')
    % plot a linear fit to see a trend
    p = polyfit(1:length(temp(1,:)),amp_post,1);
    plot(1:length(temp(1,:)), polyval(p,1:length(temp(1,:))),'k')
    xlabel('Trial #')
    ylabel({strcat('z-scored mean vessel width ',num2str(dt),'sec post-stimulus')})
    set(gca,'FontSize',14)
    box off
    
    clear t std_dev curve1 curve2 t2 temp
    
        % Here we preserve the size of the image when we save it.
    % https://dgleich.github.io/hq-matlab-figs/
    set(gcf,'InvertHardcopy','on');
    set(gcf,'PaperUnits', 'inches');
    width = 11.2; height = 8.4;
    papersize = get(gcf, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(gcf,'PaperPosition', myfiguresize);
    print(strcat('ROI',num2str(i),'_trial_width_Z'),'-dpng','-r150');
end

