function [TSdataBBBperm] = BBBpermTS(inputStacks,userInput)

inputStacks = sortedStacks;

%% create non-vascular ROI for entire x-y plane 

%go to dir w/functions
[imAn1funcDir] = getUserInput(userInput,'imAnalysis1_functions Directory');
cd(imAn1funcDir); 

%update userInput 
UIr = size(userInput,1)+1;
numROIs = input("How many BBB perm ROIs are we making? "); userInput(UIr,1) = ("How many BBB perm ROIs are we making?"); userInput(UIr,2) = (numROIs); UIr = UIr+1;

%for display purposes mostly: average across trials 
stackAVsIm = cell(1,length(inputStacks));
for z = 1:length(inputStacks)
    for trialType = 1:size(inputStacks{z},2)
        if isempty(inputStacks{z}{trialType}) == 0 
            for trial = 1:size(inputStacks{z}{trialType},2)
             stackAVsIm{z}{trialType} = mean(inputStacks{z}{trialType}{trial},3);
            end 
        end  
    end 
end 

%create the ROI boundaries           
ROIboundDatas = cell(1,numROIs);
for VROI = 1:numROIs 
    disp('Create your ROI for BBB perm analysis');

    [~,xmins,ymins,widths,heights] = firstTimeCreateROIs(1, stackAVsIm{1}{4});
    ROIboundData{1} = xmins;
    ROIboundData{2} = ymins;
    ROIboundData{3} = widths;
    ROIboundData{4} = heights;

    ROIboundDatas{VROI} = ROIboundData;
end 

%create ROIs
ROIstacks = cell(1,length(inputStacks));
for z = 1:length(inputStacks) 
    for trialType = 1:size(inputStacks{z},2)
        if isempty(inputStacks{z}{trialType}) == 0 
            for trial = 1:size(inputStacks{z}{trialType},2)
                for VROI = 1:numROIs 
                    %use the ROI boundaries to generate ROIstacks 
                    xmins = ROIboundDatas{VROI}{1};
                    ymins = ROIboundDatas{VROI}{2};
                    widths = ROIboundDatas{VROI}{3};
                    heights = ROIboundDatas{VROI}{4};
                    [ROI_stacks] = make_ROIs_notfirst_time(inputStacks{z}{trialType}{trial},xmins,ymins,widths,heights);
                    ROIstacks{z}{trialType}{trial}{VROI} = ROI_stacks;
                end        
            end 
        end 
    end 
end 

%% segment the ROIs - goal: identify non-vascular/non-terminal space 

segQ = 1; 
cd(imAn1funcDir); 
while segQ == 1     

    %segment the vessel (small sample of the data) 
    VROI = input("What BBB ROI do you want to use to make segmentation algorithm? ");

    imageSegmenter(ROIstacks{1}{4}{1}{VROI}{1}(:,:,size(ROIstacks{1}{4}{1}{VROI}{1},3)))
    continu = input('Is the image segmenter closed? Yes = 1. No = 0. ');

    while continu == 1 
        BWstacks = cell(1,length(ROIstacks));
        BW_perim = cell(1,length(ROIstacks));
        segOverlays = cell(1,length(ROIstacks));
        for Z = 1:length(ROIstacks)
            for trialType = 1:size(inputStacks{z},2)
                for VROI = 1:numROIs 
                    if isempty(inputStacks{z}{trialType}) == 0 
                        for trial = 1:size(ROIstacks{Z}{trialType},2)

                                for frame = 1:size(ROIstacks{Z}{trialType}{trial}{VROI}{1},3)
                                    [BW,~] = segmentImageBBB(ROIstacks{Z}{trialType}{trial}{VROI}{1}(:,:,frame));
                                    BWstacks{Z}{trialType}{trial}{VROI}(:,:,frame) = BW; 
                                    %get the segmentation boundaries 
                                    BW_perim{Z}{trialType}{trial}{VROI}(:,:,frame) = bwperim(BW);
                                    %overlay segmentation boundaries on data
                                    segOverlays{Z}{trialType}{trial}{VROI}(:,:,:,frame) = imoverlay(mat2gray(ROIstacks{Z}{trialType}{trial}{VROI}{1}(:,:,frame)), BW_perim{Z}{trialType}{trial}{VROI}(:,:,frame), [.3 1 .3]);
                                end 

                        end 
                    end 
                end 
            end 
        end 
        continu = 0;
    end 

    %check segmentation 
    [volIm] = getUserInput(userInput,'Is this volume imaging data? Yes = 1. Not = 0.');
    if volIm == 1
         Z = input("What Z plane do you want to see? ");
    elseif volIm == 0 
        Z = 1; 
    end 
    VROI = input("What BBB ROI do you want to see? ");
    trialType = input("What trial type do you want to see? ");
    
    %play segmentation boundaries over images 
    implay(segOverlays{Z}{trialType}{1}{VROI})

    segQ = input('Does segmentation need to be redone? Yes = 1. No = 0. ');    
end 

%% invert the mask
BWstacksInv = cell(1,length(ROIstacks));
for Z = 1:length(ROIstacks)
    for trialType = 1:size(inputStacks{z},2)
        for VROI = 1:numROIs 
            if isempty(inputStacks{z}{trialType}) == 0 
                for trial = 1:size(ROIstacks{Z}{trialType},2)
                        for frame = 1:size(ROIstacks{Z}{trialType}{trial}{VROI}{1},3)                            
                            BWstacksInv{Z}{trialType}{trial}{VROI}(:,:,frame) = ~(BWstacks{Z}{trialType}{trial}{VROI}(:,:,frame)); 
                        end 
                end 
            end 
        end 
    end 
end  

%% get pixel intensity value of extravascular space 
meanPixIntArray = cell(1,length(ROIstacks));
for Z = 1:length(ROIstacks)
    for trialType = 1:size(inputStacks{Z},2)
        for VROI = 1:numROIs 
            if isempty(inputStacks{Z}{trialType}) == 0 
                for trial = 1:size(ROIstacks{Z}{trialType},2)
                        for frame = 1:size(ROIstacks{Z}{trialType}{trial}{VROI}{1},3)                            
                            stats = regionprops(BWstacksInv{Z}{trialType}{trial}{VROI}(:,:,frame),ROIstacks{Z}{trialType}{trial}{VROI}{1}(:,:,frame),'MeanIntensity');
                            ROIpixInts = zeros(1,length(stats));
                            for stat = 1:length(stats)
                                ROIpixInts(stat) = stats(stat).MeanIntensity;
                            end 
                            meanPixIntArray{Z}{trialType}{trial}{VROI}(frame) = mean(ROIpixInts);    
                        end 
                        % turn all rows full of zeros into NaNs 
                        allZeroRows = find(all(meanPixIntArray{Z}{trialType}{trial}{VROI} == 0,2));
                        for row = 1:length(allZeroRows)
                            meanPixIntArray{Z}{trialType}{trial}{VROI} = NaN; 
                        end 
                end 
            end 
        end 
    end 
end  

           
%% normalize and z score 
            
dataMeds = cell(1,length(ROIstacks));
DFOF = cell(1,length(ROIstacks));
dataSlidingBLs = cell(1,length(ROIstacks));
DFOFsubSBLs = cell(1,length(ROIstacks));
zData = cell(1,length(ROIstacks));
for Z = 1:length(ROIstacks)
    for trialType = 1:size(inputStacks{z},2)
        for VROI = 1:numROIs 
            if isempty(inputStacks{z}{trialType}) == 0 
                for trial = 1:size(ROIstacks{Z}{trialType},2)
                    %get median value per trace
                    dataMed = nanmedian(meanPixIntArray{Z}{trialType}{trial}{VROI});     
                    dataMeds{Z}{trialType}{trial}{VROI} = dataMed;
                    %compute DF/F using median  
                    DFOF{Z}{trialType}{trial}{VROI} = (meanPixIntArray{Z}{trialType}{trial}{VROI}-dataMeds{Z}{trialType}{trial}{VROI})./dataMeds{Z}{trialType}{trial}{VROI};   
                    %get sliding baseline 
                    [dataSlidingBL]=slidingBaseline(DFOF{Z}{trialType}{trial}{VROI},floor((FPS/numZplanes)*10),0.5); %0.5 quantile thresh = the median value                 
                    dataSlidingBLs{Z}{trialType}{trial}{VROI} = dataSlidingBL;    
                    %subtract sliding baseline from DF/F
                    DFOFsubSBLs{Z}{trialType}{trial}{VROI} = DFOF{Z}{trialType}{trial}{VROI}-dataSlidingBLs{Z}{trialType}{trial}{VROI};                   
                    %z-score data                
                    zData{Z}{trialType}{trial}{VROI} = zscore(DFOFsubSBLs{Z}{trialType}{trial}{VROI});
                end 
            end 
        end 
    end 
end 

BdataToPlot = zData;

%% plot your data 
%smoothAndPlotBBBData(BdataToPlot,userInput,FPS,numZplanes,sec_before_stim_start)

%% plot averaged data 
smoothAndPlotAVBBBData(BdataToPlot,userInput,FPS,numZplanes,sec_before_stim_start)
end 