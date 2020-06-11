function [sortedData_width_Z,sortedData_width,sortedData_area,userInput,UIr,ROIboundDatas] = segmentVessels(reg_Stacks,volIm,UIr,userInput,state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,numZplanes,HDFchart)
%% create vessel segmentation ROIs - rotate if needed 
[imAn1funcDir] = getUserInput(userInput,'imAnalysis1_functions Directory');
cd(imAn1funcDir); 
numROIs = input("How many ROIs are we making? "); userInput(UIr,1) = ("How many ROIs are we making?"); userInput(UIr,2) = (numROIs); UIr = UIr+1;

rotStackAngles = zeros(1,numROIs);
ROIboundDatas = cell(1,numROIs);
for VROI = 1:numROIs 
    %rotate all the planes in Z per vessel ROI 
    [rotStacks,rotateImAngle] = rotateStack(reg_Stacks);       
    rotStackAngles(VROI) = rotateImAngle;    

    %create your ROI and apply it to all planes in Z 
    disp('Create your ROI for vessel segmentation');

    for stack = 1:length(rotStacks)   
        if stack == 1
            [ROI_stacks,xmins,ymins,widths,heights] = firstTimeCreateROIs(1,rotStacks{stack});
            ROIboundData{1} = xmins;
            ROIboundData{2} = ymins;
            ROIboundData{3} = widths;
            ROIboundData{4} = heights;
            ROIstacks{stack}{VROI} = ROI_stacks;

        elseif stack > 1 
            xmins = ROIboundData{1};
            ymins = ROIboundData{2};
            widths = ROIboundData{3};
            heights = ROIboundData{4};
            [ROI_stacks] = make_ROIs_notfirst_time(rotStacks{stack},xmins,ymins,widths,heights);
            ROIstacks{stack}{VROI} = ROI_stacks;
        end 
    end 
    ROIboundDatas{VROI} = ROIboundData;
end 


rotStackAngles = string(rotStackAngles);
rotStackAnglesJoined = join(rotStackAngles);
userInput(UIr,1) = ("ROI Rotation Angles"); userInput(UIr,2) = (rotStackAnglesJoined); UIr = UIr+1;


%% segment the vessels, get vessel width, and check segmentation
segmentVessel = 1;
BWstacks = cell(1,length(rotStacks));
BW_perim = cell(1,length(rotStacks));
segOverlays = cell(1,length(rotStacks));
frame_cnt = 500; % # of frames to use as an example in videos and segmentation validation
for VROI = 1:numROIs
    if VROI <= numROIs
        segmentVessel = 1;
    end
    while segmentVessel == 1
        %display last image of first ROI z stacks to pick the one that is most dim
        %for making segmentation algorithm
        disp('Pick a Z-stack to use for segmentation algorithm based on these pics.');
        for Zstack = 1:length(rotStacks)
            figure;
            imshow(rotStacks{Zstack}(:,:,end),[0 1800])
        end
        segm_test = input('Run Segmentation App? Yes = 1. No = 0. ');
        if segm_test==1
            dispIm = input('What Z-stack is most dim and should be used for segmentation algorithm? '); userInput(UIr,1) = ("What Z-stack is most dim and should be used for segmentation algorithm?"); userInput(UIr,2) = (dispIm); UIr = UIr+1;
            dispFr = input(strcat('What frame do you want to use for segmentation? (out of ',...
                num2str(size(ROIstacks{dispIm}{VROI}{1},3)),' frames)'));
            %segment the vessel (small sample of the data)
            imageSegmenter(ROIstacks{dispIm}{VROI}{1}(:,:,dispFr));
            continu = input('Is the image segmenter closed? Yes = 1. No = 0. ');
            segm_test = 0;
            segm_funct = uigetfile('segm*.m','Which Segmentation function do you want to use?');
        else
            continu = 1;
            segm_funct = uigetfile('segm*.m','Which Segmentation function do you want to use?');
        end
        
        while continu == 1
            % Making sure the right segmentation function was chosen - to
            % avoid crashing
            error_test = 1;
            while error_test == 1
                try
                    [test,~] = eval(strcat(segm_funct(1:end-2),'(ROIstacks{1}{1}{1}(:,:,1))'));
                    error_test = 0;
                catch
                    segm_funct = uigetfile('segm*.m','Bad choice of segmenter! Which Segmentation function do you want to use?');
                end
            end
            for Zstack = 1:length(rotStacks)
                temp = double(HDFchart(2:length(HDFchart),[1,8,9])); % Stim End Frame
                temp = temp(temp(:,3)<frame_cnt,:);
                for frame = 1:frame_cnt
                    [BW,~] = eval(strcat(segm_funct(1:end-2),'(ROIstacks{Zstack}{VROI}{1}(:,:,frame))'));
                    BWstacks{Zstack}{VROI}(:,:,frame) = BW;
                    %BWstacks{Z}{trialType}{trial}{VROI}(:,:,frame) = BW;
                    %get the segmentation boundaries
                    BW_perim{Zstack}{VROI}(:,:,frame) = bwperim(BW);
                    %overlay segmentation boundaries on data
                    if any(frame >= floor(temp(:,2)) & frame <= ceil(temp(:,3)))
                        segOverlays{Zstack}{VROI}(:,:,:,frame) = imoverlay(mat2gray(ROIstacks{Zstack}{VROI}{1}(:,:,frame)), BW_perim{Zstack}{VROI}(:,:,frame), [0.8500, 0.3250, 0.0980]);
                    else
                        segOverlays{Zstack}{VROI}(:,:,:,frame) = imoverlay(mat2gray(ROIstacks{Zstack}{VROI}{1}(:,:,frame)), BW_perim{Zstack}{VROI}(:,:,frame), [.3 1 .3]);
                    end
                end
            end
            continu = 0;
        end
        
        %check segmentation
        if volIm == 1
            Z = input("What Z plane do you want to see? ");
        elseif volIm == 0
            Z = 1;
        end
        %         VROI = input("What vessel ROI do you want to see? ");
        
        %play segmentation boundaries over images
        disp(strcat('Playing ROI',num2str(VROI)))
        implay(segOverlays{Z}{VROI})
        
        segmentVessel = input("Does the vessel need to be segmented again? Yes = 1. No = 0. ");
        if segmentVessel == 1
            clear BWthreshold BWopenRadius BW se boundaries
        end
    end
    
        %segment the vessel (all the data)
        disp('Vessel Segmentation')
        for Zstack = 1:length(rotStacks)
            %         for VROI = 1:numROIs
            for frame = 1:size(ROIstacks{1}{1}{1},3)
                [BW,~] = eval(strcat(segm_funct(1:end-2),'(ROIstacks{Zstack}{VROI}{1}(:,:,frame))'));
%                 [BW,~] = segmentImage(ROIstacks{Zstack}{VROI}{1}(:,:,frame));
                BWstacks{Zstack}{VROI}(:,:,frame) = BW;
            end
            %         end
        end
end

%% Vessel width
%get vessel width 
for Zstack = 1:length(BWstacks)
    for VROI = 1:numROIs
        [vesselDiam,bounds,area] = findVesWidth(BWstacks{Zstack}{VROI});
        area_full{Zstack}{VROI} = area; % Dmitrijs
        boundaries{Zstack}{VROI} = bounds;
        vessel_diam{Zstack}{VROI} = vesselDiam;
    end 
end 

%remove outliers = max and min value 
maxVDval = zeros(length(BWstacks),numROIs);
minVDval = zeros(length(BWstacks),numROIs);
for Zstack = 1:length(BWstacks)
    for VROI = 1:numROIs
        maxVDval(Zstack,VROI) = max(vessel_diam{Zstack}{VROI});
        minVDval(Zstack,VROI) = min(vessel_diam{Zstack}{VROI});
    end 
end 

%interpolate (average) data at frames that are max or min 
for Zstack = 1:length(BWstacks)
    for VROI = 1:numROIs 
        for frame = 2:size(ROIstacks{1}{1}{1},3)-1 %to avoid edges since i'm interpolating 
            if vessel_diam{Zstack}{VROI}(:,frame) == maxVDval(Zstack,VROI) || vessel_diam{Zstack}{VROI}(:,frame) == minVDval(Zstack,VROI)
                vessel_diam2{Zstack}{VROI}(:,frame) = (vessel_diam{Zstack}{VROI}(:,frame-1) + vessel_diam{Zstack}{VROI}(:,frame+1)) / 2 ; 
            elseif vessel_diam{Zstack}{VROI}(:,frame) ~= maxVDval(Zstack,VROI) || vessel_diam{Zstack}{VROI}(:,frame) ~= minVDval(Zstack,VROI)
                vessel_diam2{Zstack}{VROI}(:,frame) = vessel_diam{Zstack}{VROI}(:,frame);
            end 
        end 
    end 
end 

%% Creating videos and GIFs
for Zstack = 1:length(BWstacks)
    for VROI = 1:numROIs
        figure
        v = VideoWriter(strcat('segmentation_ROI',num2str(VROI),...
            '_4xSpeed_ie_',num2str(ceil(4*FPS)),'fps_',num2str(frame_cnt),'frames.mp4'),'MPEG-4');
        v.FrameRate = ceil(4*FPS);
        v.Quality = 100;
        open(v)
        gif_file = strcat('segmentation_ROI',num2str(VROI),...
            '_4xSpeed_ie_',num2str(ceil(4*FPS)),'fps_',num2str(frame_cnt),'frames.gif');
        for frame = 1:500
            imagesc(segOverlays{Zstack}{VROI}(:,:,:,frame));
            % Finding approx middle of bounding box
            temp = size(segOverlays{Zstack}{VROI}(:,:,:,frame),1)/2;
            x = [find(BW_perim{Zstack}{VROI}(ceil(temp),:,frame),1,'first') 
                find(BW_perim{Zstack}{VROI}(ceil(temp),:,frame),1,'last')];
            x = mean(x)-vessel_diam2{Zstack}{VROI}(frame)/2 :1 : ...
                mean(x)+vessel_diam2{Zstack}{VROI}(frame)/2;
            y = zeros(1,length(x));
            y(1:end) = size(segOverlays{Zstack}{VROI}(:,:,:,frame),1)/2;
            if (length(x)>1) % to avoid the error during 1st frame when there is no measurement
                line(x,y,'LineWidth',8,'Color','red')
                text(mean(x)-0.5,y(1)-0.5,num2str(round(vessel_diam2{Zstack}{VROI}(frame),1)),'FontSize',10,'Color','r')
            end
            F = getframe(gcf);
            writeVideo(v,F)
            % Write to the GIF File
            [A,map] = rgb2ind(frame2im(F),256);
            if frame == 1
                imwrite(A,map,gif_file,'gif', 'Loopcount',inf,'DelayTime',1/(FPS*4));
            else
                imwrite(A,map,gif_file,'gif','WriteMode','append','DelayTime',1/(FPS*4));
            end
        end
        close(v)
    end
end

%% Normalized diameter estimation
dataMeds = cell(1,length(BWstacks));
DFOF = cell(1,length(BWstacks));
dataSlidingBLs = cell(1,length(BWstacks));
DFOFsubSBLs = cell(1,length(BWstacks));
zVData = cell(1,length(BWstacks));
vessel_diam3 = cell(1,length(BWstacks));
area_full2 = cell(1,length(BWstacks));
for z = 1:length(BWstacks)
    for VROI = 1:numROIs  
        vessel_diam3{z}(VROI,:) = vessel_diam2{z}{VROI};
        area_full2{z}(VROI,:) = area_full{z}{VROI};
        %get median value per trace
        dataMed = median(vessel_diam2{z}{VROI});     
        dataMeds{z}(VROI,:) = dataMed;        
        %compute DF/F using means  
        DFOF{z}(VROI,:) = (vessel_diam2{z}{VROI}-dataMeds{z}(VROI,:))./dataMeds{z}(VROI,:);              
        %get sliding baseline 
        [dataSlidingBL]=slidingBaseline(DFOF{z}(VROI,:),floor((FPS/numZplanes)*10),0.5); %0.5 quantile thresh = the median value                 
        dataSlidingBLs{z}(VROI,:) = dataSlidingBL;     
        %subtract sliding baseline from DF/F
        DFOFsubSBLs{z}(VROI,:) = DFOF{z}(VROI,:)-dataSlidingBLs{z}(VROI,:);        
        %z-score data 
        zVData{z}(VROI,:) = zscore(DFOFsubSBLs{z}(VROI,:));
    end
end


%% separate data by trial type 
disp('Organizing Data by Trial Type')

%separate the data 
sortedData_width_Z = cell(1,length(BWstacks));
sortedData_width = cell(1,length(BWstacks));
sortedData_area = cell(1,length(BWstacks));
for Zstack = 1:length(BWstacks)
  for VROI = 1:numROIs
      [sorted_Data] = eventTriggeredAverages(zVData{Zstack}(VROI,:),state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
      sortedData_width_Z{Zstack}{VROI} = sorted_Data;
      clear sorted_Data
      [sorted_Data] = eventTriggeredAverages(vessel_diam3{Zstack}(VROI,:),state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
      sortedData_width{Zstack}{VROI} = sorted_Data;
      clear sorted_Data
      [sorted_Data] = eventTriggeredAverages(area_full2{Zstack}(VROI,:),state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes);
      sortedData_area{Zstack}{VROI} = sorted_Data;
      clear sorted_Data
  end 
end 

end 
