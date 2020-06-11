function [regStacks,userInput,UIr,state_start_f,state_end_f,vel_wheel_data,TrialTypes,HDFchart] = imRegistration(imAnDirs)
UIr = 1;
userInput(UIr,1) = ("imAnalysis1_functions Directory"); userInput(UIr,2) = (imAnDirs(1,:)); UIr = UIr+1;
userInput(UIr,1) = ("imAnalysis2_functions Directory"); userInput(UIr,2) = (imAnDirs(2,:)); UIr = UIr+1;

volIm = input("Is this volume imaging data? Yes = 1. No = 0. "); userInput(UIr,1) = ("Is this volume imaging data? Yes = 1. Not = 0."); userInput(UIr,2) = (volIm); UIr = UIr+1;
if volIm == 1
    splitType = input('How must the data be split? Serial Split = 0. Alternating Split = 1. '); userInput(UIr,1) = ("How must the data be split? Serial Split = 0. Alternating Split = 1."); userInput(UIr,2) = (splitType); UIr = UIr+1;
    numZplanes = input('How many Z planes are there? '); userInput(UIr,1) = ("How many Z planes are there?"); userInput(UIr,2) = (numZplanes); UIr = UIr+1;
elseif volIm == 0
    numZplanes = 1; userInput(UIr,1) = ("How many Z planes are there?"); userInput(UIr,2) = (numZplanes); UIr = UIr+1;
end 
framePeriod = input("What is the framePeriod? "); userInput(UIr,1) = ("What is the framePeriod? "); userInput(UIr,2) = (framePeriod); UIr = UIr+1;
state = input("What teensy state does the stimulus happen in? "); userInput(UIr,1) = ("What teensy state does the stimulus happen in?"); userInput(UIr,2) = (state); UIr = UIr+1;
stimNumQ = input("Is there more than one stimulus type? Yes = 1. No = 0."); userInput(UIr,1) = ("Is there more than one stimulus type? Yes = 1. No = 0."); userInput(UIr,2) = (stimNumQ); UIr = UIr+1; 
if stimNumQ == 1
    stimNums = input("How many different kinds of stimuli were used?"); userInput(UIr,1) = ("How many different kinds of stimuli were used?"); userInput(UIr,2) = (stimNums); UIr = UIr+1; 
    stimTypes = cell(1,stimNums);
    for stim = 1:stimNums
        stimTypes{stim} = input(sprintf("What type of stimulus is Stim #%d? (put strings in quotes) ",stim)); 
        userInput(UIr,1) = (sprintf("Stim Type #%d",stim)); userInput(UIr,2) = (stimTypes{stim}); UIr = UIr+1;     
    end  
end 
stimLengthQ = input("Are all stim trains the same length in time? Yes = 1. No = 0. "); userInput(UIr,1) = ("Are all stim trains the same length in time?"); userInput(UIr,2) = (stimLengthQ); UIr = UIr+1; 
if stimLengthQ == 1
    stimLengthsQ = 1; userInput(UIr,1) = ("How many different stimulus time lengths are there?"); userInput(UIr,2) = (stimLengthsQ); UIr = UIr+1; 
    stimTimes = zeros(1,stimLengthsQ);
    for stimTime = 1:stimLengthsQ
        stimTimes(stimTime) = input(sprintf("How long is stim length #%d in seconds? ",stimTime)); 
    end 
    stimTimes = string(stimTimes);
    stimTimesJoined = join(stimTimes);
    userInput(UIr,1) = ("Stim Time Lengths (sec)"); userInput(UIr,2) = (stimTimesJoined); UIr = UIr+1;     
elseif stimLengthQ == 0
    stimLengthsQ = input("How many different stimulus time lengths are there? "); userInput(UIr,1) = ("How many different stimulus time lengths are there?"); userInput(UIr,2) = (stimLengthsQ); UIr = UIr+1; 
    stimTimes = zeros(1,stimLengthsQ);    
    for stimTime = 1:stimLengthsQ
        stimTimes(stimTime) = input(sprintf("How long is stim length #%d in seconds? ",stimTime)); 
    end 
    stimTimes = string(stimTimes);
    stimTimesJoined = join(stimTimes);
    userInput(UIr,1) = ("Stim Time Lengths (sec)"); userInput(UIr,2) = (stimTimesJoined); UIr = UIr+1;       
end 
channelOfInterest = input('What color channel is most relevant? Red = 0. Green = 1. '); userInput(UIr,1) = ("What color channel is most relevant?"); userInput(UIr,2) = (channelOfInterest); UIr = UIr+1;
%% make HDF chart and get wheel data  
disp('Making HDF Chart')
[imAn1funcDir] = getUserInput(userInput,'imAnalysis1_functions Directory');
cd(imAn1funcDir);
[HDFchart,state_start_f,state_end_f,FPS,vel_wheel_data,TrialTypes] = makeHDFchart_redBlueStim(state,framePeriod);
userInput(UIr,1) = ("FPS"); userInput(UIr,2) = (FPS); UIr = UIr+1; 

%% get the images 
disp('Importing Images')
cd(uigetdir('*.*','WHERE ARE THE PHOTOS? '));    
fileList = dir('*.tif*');
imStackDir = fileList(1).folder; %userInput(UIr,1) = ("Where Photos Are Found"); userInput(UIr,2) = (imStackDir); UIr = UIr+1;
greenStackLength = size(fileList,1)/2;
redStackLength = size(fileList,1)/2;
image = imread(fileList(1).name); 
redImageStack = zeros(size(image,1),size(image,2),redStackLength);
for frame = 1:redStackLength 
    image = imread(fileList(frame).name); 
    redImageStack(:,:,frame) = image;
end 
count = 1; 
greenImageStack = zeros(size(image,1),size(image,2),greenStackLength);
for frame = redStackLength+1:(greenStackLength*2) 
    image = imread(fileList(frame).name); 
    greenImageStack(:,:,count) = image;
    count = count + 1;
end 

%% separate volume imaging data into separate stacks per z plane and motion correction
if volIm == 1
    disp('Separating Z-planes')
    %reorganize data by zPlane and prep for motion correction 
    [gZstacks] = splitStacks(greenImageStack,splitType,numZplanes);
    [gVolStack] = reorgVolStack(greenImageStack,splitType,numZplanes);
    [rZstacks] = splitStacks(redImageStack,splitType,numZplanes);
    [rVolStack] = reorgVolStack(redImageStack,splitType,numZplanes);    
    
    if channelOfInterest == 1
        %2D rigid registration
        disp('2D Motion Correction')
        gRegTemplates = cell(1,size(gZstacks,2));
        ggRegZstacks = cell(1,size(gZstacks,2));
%         rRegTemplates = cell(1,size(gZstacks,2));
%         grRegZstacks = cell(1,size(gZstacks,2));
        for Zstack = 1:size(gZstacks,2)    
            gRegTemplates{Zstack} = mean(gZstacks{Zstack},3);
            [ggRegStack,~] = registerVesStack(gZstacks{Zstack},gRegTemplates{Zstack}); % the seconds output = transformations - I'm not currently using these 
            ggRegZstacks{Zstack} = ggRegStack;

            rRegTemplates{Zstack} = mean(rZstacks{Zstack},3);
            [grRegStack,~] = registerVesStack(gZstacks{Zstack},rRegTemplates{Zstack}); % the seconds output = transformations - I'm not currently using these 
            grRegZstacks{Zstack} = grRegStack;
        end              

        %3D registration     
        disp('3D Motion Correction')
        %need minimum 4 planes in Z for 3D registration to work-time to interpolate
        gVolStack5 = zeros(size(gVolStack,1),size(gVolStack,2),size(gVolStack,3)+size(gVolStack,3)-1,size(gVolStack,4));
%         rVolStack5 = zeros(size(rVolStack,1),size(rVolStack,2),size(gVolStack,3)+size(gVolStack,3)-1,size(gVolStack,4));
        for ind = 1:size(gVolStack,4)
            count = 1;
            for zplane = 1:size(gVolStack,3)+size(gVolStack,3)-1
                if rem(zplane,2) == 0
                    gVolStack5(:,:,zplane,ind) = mean(gVolStack(:,:,count-1:count,ind),3);
%                     rVolStack5(:,:,zplane,ind) = mean(rVolStack(:,:,count-1:count,ind),3);
                elseif rem(zplane,2) == 1 
                    gVolStack5(:,:,zplane,ind) = gVolStack(:,:,count,ind);
%                     rVolStack5(:,:,zplane,ind) = rVolStack(:,:,count,ind);
                    count = count +1;
                end 
            end 
        end 

        gTemplate = mean(gVolStack5,4);
%         rTemplate = mean(rVolStack5,4);
        %create optimizer and metric, setting modality to multimodal because the
        %template and actual images look different (as template is mean = smoothed)
        [optimizer, metric] = imregconfig('multimodal')    ;
        %tune the poperties of the optimizer
        optimizer.InitialRadius = 0.004;
        optimizer.Epsilon = 1.5;
        optimizer.GrowthFactor = 1.01;
        optimizer.MaximumIterations = 300;
        for ind = 1:size(gVolStack,4)
            ggRegVolStack(:,:,:,ind) = imregister(gVolStack5(:,:,:,ind),gTemplate,'affine', optimizer, metric,'PyramidLevels',1);
%             grRegVolStack(:,:,:,ind) = imregister(gVolStack5(:,:,:,ind),rTemplate,'affine', optimizer, metric,'PyramidLevels',1);
        end 
        [ggVZstacks5] = volStack2splitStacks(ggRegVolStack);
%         [grVZstacks5] = volStack2splitStacks(grRegVolStack);
   
        %get rid of extra z planes
        count = 1; 
        ggVZstacks3 = cell(1,size(gZstacks,2));
%         grVZstacks3 = cell(1,size(gZstacks,2));
        for zPlane = 1:size(gZstacks,2)
            ggVZstacks3{zPlane} = ggVZstacks5{count};
%             grVZstacks3{zPlane} = grVZstacks5{count};
            count = count+2;
        end 
        
        %package data for output 
%         regStacks{2,1} = ggRegZstacks; regStacks{1,1} = 'ggRegZstacks';
%         regStacks{2,2} = grRegZstacks; regStacks{1,2} = 'grRegZstacks';
        regStacks{2,3} = ggVZstacks3; regStacks{1,3} = 'ggVZstacks3';
%         regStacks{2,4} = grVZstacks3; regStacks{1,4} = 'grVZstacks3';
        
           
    elseif channelOfInterest == 0 
%         %2D rigid registration
%         disp('2D Motion Correction')
%         gRegTemplates = cell(1,size(rZstacks,2));
%         rRegTemplates = cell(1,size(rZstacks,2));
%         rgRegZstacks = cell(1,size(rZstacks,2));
%         rrRegZstacks = cell(1,size(rZstacks,2));
%         for Zstack = 1:size(rZstacks,2)    
%             gRegTemplates{Zstack} = mean(gZstacks{Zstack},3);
%             [rgRegStack,~] = registerVesStack(rZstacks{Zstack},gRegTemplates{Zstack}); % the seconds output = transformations - I'm not currently using these 
%             rgRegZstacks{Zstack} = rgRegStack;
% 
%             rRegTemplates{Zstack} = mean(rZstacks{Zstack},3);
%             [rrRegStack,~] = registerVesStack(rZstacks{Zstack},rRegTemplates{Zstack}); % the seconds output = transformations - I'm not currently using these 
%             rrRegZstacks{Zstack} = rrRegStack;
%         end              
        
        %3D registration     
        disp('3D Motion Correction')
        %need minimum 4 planes in Z for 3D registration to work-time to interpolate
%         gVolStack5 = zeros(size(gVolStack,1),size(gVolStack,2),size(gVolStack,3)+size(gVolStack,3)-1,size(gVolStack,4));
        rVolStack5 = zeros(size(rVolStack,1),size(rVolStack,2),size(gVolStack,3)+size(gVolStack,3)-1,size(gVolStack,4));
        for ind = 1:size(gVolStack,4)
            count = 1;
            for zplane = 1:size(gVolStack,3)+size(gVolStack,3)-1
                if rem(zplane,2) == 0
%                     gVolStack5(:,:,zplane,ind) = mean(gVolStack(:,:,count-1:count,ind),3);
                    rVolStack5(:,:,zplane,ind) = mean(rVolStack(:,:,count-1:count,ind),3);
                elseif rem(zplane,2) == 1 
%                     gVolStack5(:,:,zplane,ind) = gVolStack(:,:,count,ind);
                    rVolStack5(:,:,zplane,ind) = rVolStack(:,:,count,ind);
                    count = count +1;
                end 
            end 
        end 

%         gTemplate = mean(gVolStack5,4);
        rTemplate = mean(rVolStack5,4);
        %create optimizer and metric, setting modality to multimodal because the
        %template and actual images look different (as template is mean = smoothed)
        [optimizer, metric] = imregconfig('multimodal')    ;
        %tune the poperties of the optimizer
        optimizer.InitialRadius = 0.004;
        optimizer.Epsilon = 1.5;
        optimizer.GrowthFactor = 1.01;
        optimizer.MaximumIterations = 300;
        parfor ind = 1:size(gVolStack,4)
%             rgRegVolStack(:,:,:,ind) = imregister(rVolStack5(:,:,:,ind),gTemplate,'affine', optimizer, metric,'PyramidLevels',1);
            rrRegVolStack(:,:,:,ind) = imregister(rVolStack5(:,:,:,ind),rTemplate,'affine', optimizer, metric,'PyramidLevels',1);
        end 
%         [rgVZstacks5] = volStack2splitStacks(rgRegVolStack);
        [rrVZstacks5] = volStack2splitStacks(rrRegVolStack);
   
        %get rid of extra z planes
        count = 1; 
%         rgVZstacks3 = cell(1,size(gZstacks,2));
        rrVZstacks3 = cell(1,size(gZstacks,2));
        for zPlane = 1:size(gZstacks,2)
%             rgVZstacks3{zPlane} = rgVZstacks5{count};
            rrVZstacks3{zPlane} = rrVZstacks5{count};
            count = count+2;
        end 
        
        %package data for output 
%         regStacks{2,1} = rgRegZstacks; regStacks{1,1} = 'rgRegZstacks';
%         regStacks{2,2} = rrRegZstacks; regStacks{1,2} = 'rrRegZstacks';
%         regStacks{2,3} = rgVZstacks3; regStacks{1,3} = 'rgVZstacks3';
        regStacks{2,4} = rrVZstacks3; regStacks{1,4} = 'rrVZstacks3';
    end 
            
elseif volIm == 0
    if channelOfInterest == 1
        disp('2D Motion Correction')
        %2D register imaging data    
        gTemplate = mean(greenImageStack,3);
        [ggRegStack,~] = registerVesStack(greenImageStack,gTemplate);  
        ggRegZstacks{1} = ggRegStack;
%         rTemplate = mean(redImageStack,3);
%         [grRegStack,~] = registerVesStack(greenImageStack,rTemplate);  
%         grRegZstacks{1} = grRegStack;
%         
        %package data for output 
        regStacks{2,1} = ggRegZstacks; regStacks{1,1} = 'ggRegZstacks';
%         regStacks{2,2} = grRegZstacks; regStacks{1,2} = 'grRegZstacks';
   
    elseif channelOfInterest == 0
        disp('2D Motion Correction')
        %2D register imaging data    
%         gTemplate = mean(greenImageStack,3);
%         [rgRegStack,~] = registerVesStack(redImageStack,gTemplate);  
%         rgRegZstacks{1} = rgRegStack;
        rTemplate = mean(redImageStack,3);
        [rrRegStack,~] = registerVesStack(redImageStack,rTemplate);  
        rrRegZstacks{1} = rrRegStack;
        
        %package data for output 
%         regStacks{2,1} = rgRegZstacks; regStacks{1,1} = 'rgRegZstacks';
        regStacks{2,2} = rrRegZstacks; regStacks{1,2} = 'rrRegZstacks';
    end     
end 

%% check registration 

if volIm == 1 
    if channelOfInterest == 1
        %check relationship b/w template and 2D registered images
        ggTemp2regCorr2D = zeros(size(gZstacks,2),size(ggVZstacks3{1},3));
        grTemp2regCorr2D = zeros(size(gZstacks,2),size(ggVZstacks3{1},3));
        for zPlane = 1:size(gZstacks,2)
            for ind = 1:size(ggVZstacks3{1},3)
                ggTemp2regCorr2D(zPlane,ind) = corr2(gRegTemplates{zPlane},ggRegZstacks{zPlane}(:,:,ind));
                grTemp2regCorr2D(zPlane,ind) = corr2(rRegTemplates{zPlane},grRegZstacks{zPlane}(:,:,ind));
            end 
        end 
        
        %plot 2D registrations for comparison 
        figure;
        subplot(1,2,1);
        hold all; 
        for zPlane = 1:size(gZstacks,2)
            plot(ggTemp2regCorr2D(zPlane,:));
            title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Green Channel Registered with Green Channel Template'}); 
        end 
        subplot(1,2,2);
        hold all; 
        for zPlane = 1:size(gZstacks,2)
            plot(grTemp2regCorr2D(zPlane,:));
            title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Green Channel Registered with Red Channel Template'}); 
        end 
    
        
        %check relationship b/w template and 3D registered images
        count = 1; 
        ggTemp2regCorr3D = zeros(size(gZstacks,2),size(ggVZstacks3{1},3));
%         grTemp2regCorr3D = zeros(size(gZstacks,2),size(ggVZstacks3{1},3));
        for zPlane = 1:size(gZstacks,2)
            for ind = 1:size(ggVZstacks3{1},3)
                ggTemp2regCorr3D(zPlane,ind) = corr2(gTemplate(:,:,count),ggVZstacks3{zPlane}(:,:,ind));
%                 grTemp2regCorr3D(zPlane,ind) = corr2(rTemplate(:,:,count),grVZstacks3{zPlane}(:,:,ind));
            end 
            count = count+2;
        end 
        
        %plot 3D registration for comparison 
        figure;
%         subplot(1,2,1);
        hold all; 
        for zPlane = 1:size(gZstacks,2)
            plot(ggTemp2regCorr3D(zPlane,:));
            title({'Correlation Coefficient of 3D Motion Correction Template and Output';'Green Channel Registered with Green Channel Template'}); 
        end 
%         subplot(1,2,2);
%         hold all; 
%         for zPlane = 1:size(gZstacks,2)
%             plot(grTemp2regCorr3D(zPlane,:));
%             title({'Correlation Coefficient of 3D Motion Correction Template and Output';'Green Channel Registered with Red Channel Template'}); 
%         end 
    
   elseif channelOfInterest == 0
%         %         %check relationship b/w template and 2D registered images
% %         rgTemp2regCorr2D = zeros(size(gZstacks,2),size(rgVZstacks3{1},3));
%         rrTemp2regCorr2D = zeros(size(gZstacks,2),size(rrVZstacks3{1},3));
%         for zPlane = 1:size(gZstacks,2)
%             for ind = 1:size(rrVZstacks3{1},3)
% %                 rgTemp2regCorr2D(zPlane,ind) = corr2(gRegTemplates{zPlane},rgRegZstacks{zPlane}(:,:,ind));
%                 rrTemp2regCorr2D(zPlane,ind) = corr2(rRegTemplates{zPlane},rrRegZstacks{zPlane}(:,:,ind));
%             end 
%         end 

%         %plot 2D registrations for comparison 
%         figure;
%         subplot(1,2,1);
%         hold all; 
%         for zPlane = 1:size(gZstacks,2)
%             plot(rgTemp2regCorr2D(zPlane,:));
%             title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Red Channel Registered with Green Channel Template'}); 
%         end 
%         subplot(1,2,2);
%         hold all; 
%         for zPlane = 1:size(gZstacks,2)
%             plot(rrTemp2regCorr2D(zPlane,:));
%             title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Red Channel Registered with Red Channel Template'}); 
%         end 
 
        %check relationship b/w template and 3D registered images
        count = 1; 
%         rgTemp2regCorr3D = zeros(size(gZstacks,2),size(rgVZstacks3{1},3));
        rrTemp2regCorr3D = zeros(size(rZstacks,2),size(rrVZstacks3{1},3));
        for zPlane = 1:size(rZstacks,2)
            for ind = 1:size(rrVZstacks3{1},3)
%                 rgTemp2regCorr3D(zPlane,ind) = corr2(gTemplate(:,:,count),rgVZstacks3{zPlane}(:,:,ind));
                rrTemp2regCorr3D(zPlane,ind) = corr2(rTemplate(:,:,count),rrVZstacks3{zPlane}(:,:,ind));
            end 
            count = count+2;
        end 
        
        %plot 3D registrations for comparison 
        figure;
%         subplot(1,2,1);
%         hold all; 
%         for zPlane = 1:size(gZstacks,2)
%             plot(rgTemp2regCorr3D(zPlane,:));
%             title({'Correlation Coefficient of 3D Motion Correction Template and Output';'Red Channel Registered with Green Channel Template'}); 
%         end 
%         subplot(1,2,2);
        hold all; 
        for zPlane = 1:size(rZstacks,2)
            plot(rrTemp2regCorr3D(zPlane,:));
            title({'Correlation Coefficient of 3D Motion Correction Template and Output';'Red Channel Registered with Red Channel Template'}); 
        end 
   end    
elseif volIm == 0 
    
    if channelOfInterest == 1
        %check relationship b/w template and 2D registered images     
        ggTemp2regCorr2D = zeros(1,size(ggRegZstacks{1},3));
%         grTemp2regCorr2D = zeros(1,size(ggRegZstacks{1},3));
        for ind = 1:size(ggRegZstacks{1},3)
            ggTemp2regCorr2D(ind) = corr2(gTemplate,ggRegZstacks{1}(:,:,ind));
%             grTemp2regCorr2D(ind) = corr2(rTemplate,grRegZstacks{1}(:,:,ind));
        end 
        
        
        %plot 2D registrations for comparison 
        figure;
%         subplot(1,2,1);
        plot(ggTemp2regCorr2D);
        title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Green Channel Registered with Green Channel Template'}); 
%         subplot(1,2,2);
% %         plot(grTemp2regCorr2D);
%         title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Green Channel Registered with Red Channel Template'}); 

    
   elseif channelOfInterest == 0
        %check relationship b/w template and 2D registered images
%         rgTemp2regCorr2D = zeros(1,size(rgRegZstacks{1},3));
        rrTemp2regCorr2D = zeros(1,size(rrRegZstacks{1},3));
        for ind = 1:size(rrRegZstacks{1},3)
%             rgTemp2regCorr2D(ind) = corr2(gTemplate,rgRegZstacks{1}(:,:,ind));
            rrTemp2regCorr2D(ind) = corr2(rTemplate,rrRegZstacks{1}(:,:,ind));
        end 

        %plot 2D registrations for comparison 
        figure;
%         subplot(1,2,1);
%         plot(rgTemp2regCorr2D);
%         title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Red Channel Registered with Green Channel Template'}); 
%         subplot(1,2,2);
        plot(rrTemp2regCorr2D);
        title({'Correlation Coefficient of 2D Motion Correction Template and Output';'Red Channel Registered with Red Channel Template'}); 
   end            
end 
end 