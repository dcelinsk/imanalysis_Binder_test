%get just the wheel data
temp = matfile('SFWT6_20190607_ROI4_BBB.mat');
SFWT6_ROI4_Wdata = temp.wheelDataToPlot;

%% resample 
[RSFWT6_ROI3_Wdata,RSFWT6_ROI4_Wdata] = resampleWheelData(SFWT6_ROI3_Wdata,SFWT6_ROI4_Wdata);
% [RSF57_ROI1_Wdata,RSF57_ROI2_Wdata] = resampleWheelData(SF57_ROI1_Wdata,SF57_ROI2_Wdata);
% [RSF56_ROI1_Wdata,RSF56_ROI2_Wdata] = resampleWheelData(SF56_ROI1_Wdata,SF56_ROI2_Wdata);
% [RSF53_ROI1_Wdata,RSF53_ROI2_Wdata] = resampleWheelData(SF53_ROI1_Wdata,SF53_ROI2_Wdata);
% RSF63_ROI1_Wdata = SF63_ROI1_Wdata;
% RSF64_ROI1_Wdata = SF64_ROI1_Wdata;

%% average across trials
WAVdataToPlot = RSFWT6_ROI4_Wdata;

%average across trials
AVarray = cell(1,length(WAVdataToPlot));
AVdata = cell(1,length(WAVdataToPlot));
for trialType = 1:size(WAVdataToPlot,2)
    if isempty(WAVdataToPlot{trialType}) == 0 
        for trial = 1:size(WAVdataToPlot{trialType},2)
            if isempty(WAVdataToPlot{trialType}{trial}) == 0 
                AVarray{trialType}(trial,:) = WAVdataToPlot{trialType}{trial};
                AVdata{trialType} = nanmean(AVarray{trialType},1);
            end 
        end 
    end   
end    

SFWT6_ROI4av = AVdata;

clear AVarray AVdata

%% average across imaging FOVs 

for trialType = 1:size(SFWT6_ROI3av,2)
    SFWT6av{trialType} = (SFWT6_ROI3av{trialType} + SFWT6_ROI4av{trialType})/2;
%     SF57av{trialType} = (SF57_ROI1av{trialType} + SF57_ROI2av{trialType})/2;
%     SF56av{trialType} = (SF56_ROI1av{trialType} + SF56_ROI2av{trialType})/2;
%     SF53av{trialType} = (SF53_ROI1av{trialType} + SF53_ROI2av{trialType})/2;
end 
% SF63av = SF63_ROI1av;
% SF64av = SF64_ROI1av;





%% resample averaged individual mouse data 
for trialType = 1:size(SF63_ROI1av,2)
    if trialType == 1 || trialType == 3 
        goalLen = length(SF64av{1});
        RSF63av{trialType} = resample(SF63av{trialType},goalLen,length(SF63av{trialType}));
%         RSF57av{trialType} = resample(SF57av{trialType},goalLen,length(SF57av{trialType}));
%         RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
    elseif trialType == 2 || trialType == 4 
        goalLen = length(SF64av{4});
        RSF63av{trialType} = resample(SF63av{trialType},goalLen,length(SF63av{trialType}));
%         RSF57av{trialType} = resample(SF57av{trialType},goalLen,length(SF57av{trialType}));
%         RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
    end 
end 

RSF64av = SF64av;
RSF61av = SF61av;

%% put all data into the same array 
% for trialType = 1:size(SF64_ROI1av,2)
% %     if isempty(SF56av{trialType}) == 0 
%         miceData{trialType}(1,:) = RSF63av{trialType};
%         miceData{trialType}(2,:) = RSF64av{trialType};
% %         miceData{trialType}(3,:) = RSF58av{trialType};
% %         miceData{trialType}(4,:) = RSF53av{trialType};
% %     end 
% end 
% 
% %% get SEM and average across mice 
% for trialType = 1:size(SF63_ROI1av,2)
%     varMiceData{trialType} = nanvar(miceData{trialType});
%     avMiceData{trialType} = nanmean(miceData{trialType});
% end 
miceData = SFWT6av;
avMiceData = SFWT6av;
%% smooth data if you want 
FPS = length(miceData{4})/60;
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ == 1 
    filtTime = input('How many seconds do you want to smooth your data by? '); 
    
    for trialType = 1:size(avMiceData,2)   
        if isempty(avMiceData{trialType}) == 0                  
            [VfiltD] = MovMeanSmoothData(avMiceData{trialType},filtTime,FPS);
%             [VfiltV] = MovMeanSmoothData(varMiceData{trialType},filtTime,FPS);
            VfiltData{trialType} = VfiltD;   
%             VfiltVar{trialType} = VfiltV; 
        end 
    end
     
elseif smoothQ == 0

    for trialType = 1:size(avMiceData,2)   
        if isempty(avMiceData{trialType}) == 0                           
            VfiltData{trialType} = avMiceData{trialType};   
%             VfiltVar{trialType} = varMiceData{trialType};
        end 
    end
end 

for trialType = 1:size(avMiceData,2)
    semMiceData{trialType} = ((nanstd(miceData{trialType},1)))/(sqrt(size(miceData{trialType},1)));
end 

%% plot 
% dataMin = input("data Y axis MIN: ");
% dataMax = input("data Y axis MAX: ");
FPSstack = FPS;
baselineEndFrame = round(20*(FPSstack));



for trialType = 1:size(avMiceData,2)  
    if isempty(avMiceData{trialType}) == 0

        figure;
        %set time in x axis            
        if trialType == 1 || trialType == 3 
            Frames = size(VfiltData{trialType},2);                
            Frames_pre_stim_start = -((Frames-1)/2); 
            Frames_post_stim_start = (Frames-1)/2; 
            sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+1);
            FrameVals = round((1:FPSstack*2:Frames)-1); 
        elseif trialType == 2 || trialType == 4 
            Frames = size(VfiltData{trialType},2);
            Frames_pre_stim_start = -((Frames-1)/2); 
            Frames_post_stim_start = (Frames-1)/2; 
            sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+10);
            FrameVals = round((1:FPSstack*2:Frames)-1); 
        end 
        ax=gca;
        plot(VfiltData{trialType},'k','LineWidth',2)
        hold all;     

%         varargout = boundedline(1:size(VfiltData{trialType},2),VfiltData{trialType},semMiceData{trialType},'k','transparency', 0.3,'alpha');                                                                             

        ax.XTick = FrameVals;
        ax.XTickLabel = sec_TimeVals;
        ax.FontSize = 20;
        if trialType == 1 
            plot([round(baselineEndFrame+((FPSstack)*2)) round(baselineEndFrame+((FPSstack)*2))], [-5000 5000], 'b','LineWidth',3)
            %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*2)) round(baselineEndFrame+((FPSstack/numZplanes)*2)) baselineEndFrame],[-5000 -5000 5000 5000],'b')
            %alpha(0.4)   
            plot([baselineEndFrame baselineEndFrame], [-5000 5000], 'b','LineWidth',3) 
        elseif trialType == 3 
            plot([round(baselineEndFrame+((FPSstack)*2)) round(baselineEndFrame+((FPSstack)*2))], [-5000 5000], 'r','LineWidth',3)
            %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*2)) round(baselineEndFrame+((FPSstack/numZplanes)*2)) baselineEndFrame],[-5000 -5000 5000 5000],'r')
            %alpha(0.4)     
            plot([baselineEndFrame baselineEndFrame], [-5000 5000], 'r','LineWidth',3) 
        elseif trialType == 2 
            plot([round(baselineEndFrame+((FPSstack)*20)) round(baselineEndFrame+((FPSstack)*20))], [-5000 5000], 'b','LineWidth',3)
            %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*20)) round(baselineEndFrame+((FPSstack/numZplanes)*20)) baselineEndFrame],[-5000 -5000 5000 5000],'b')
            %alpha(0.4)   
            plot([baselineEndFrame baselineEndFrame], [-5000 5000], 'b','LineWidth',3) 
        elseif trialType == 4 
            plot([round(baselineEndFrame+((FPSstack)*20)) round(baselineEndFrame+((FPSstack)*20))], [-5000 5000], 'r','LineWidth',3)
            %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*20)) round(baselineEndFrame+((FPSstack/numZplanes)*20)) baselineEndFrame],[-5000 -5000 5000 5000],'r')
            %alpha(0.4)  
            plot([baselineEndFrame baselineEndFrame], [-5000 5000], 'r','LineWidth',3) 
        end

%         legend('BBB','DA Calcium','Vessel Width')
%             legend('BBB','Vessel Width')
        ylim([dataMin dataMax]);

        if smoothQ == 1 
            title(sprintf('Wheel data across mice smoothed by %d seconds.',filtTime));
        elseif smoothQ == 0
            title("Raw wheel data across mice.");
        end 

    end                       
end

