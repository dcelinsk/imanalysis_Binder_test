%% get just the V data
temp = matfile('SFWT6_20190607_ROI4_Vwidth.mat');
SFWT6_ROI4_Vdata = temp.VdataToPlot;

%% resample 
% [RSF58_ROI1_Vdata,RSF58_ROI2_Vdata] = resampleVWdata(SF58_ROI1_Vdata,SF58_ROI2_Vdata);
[RSFWT6_ROI3_Vdata,RSFWT6_ROI4_Vdata] = resampleVWdata(SFWT6_ROI3_Vdata,SFWT6_ROI4_Vdata);
% [RSF56_ROI1_Vdata,RSF56_ROI2_Vdata] = resampleVWdata(SF56_ROI1_Vdata,SF56_ROI2_Vdata);
% [RSF53_ROI1_Vdata,RSF53_ROI2_Vdata] = resampleVWdata(SF53_ROI1_Vdata,SF53_ROI2_Vdata);
% RSF64_ROI1_Vdata = SF64_ROI1_Vdata;

%% average across planes in Z, ROIs, and trials  
VdataToPlot = RSFWT6_ROI4_Vdata;

for z = 1:length(VdataToPlot)
    for ROI = 1:size(VdataToPlot{z},2)
        for trialType = 1:size(VdataToPlot{1}{1},2)   
            if isempty(VdataToPlot{z}{ROI}{trialType}) == 0                  
                for trial = 1:length(VdataToPlot{z}{ROI}{trialType})      
                     if isempty(VdataToPlot{z}{ROI}{trialType}{trial}) == 0    
                        
                        VAVdataToPlot1_array{z}{trialType}{trial}(ROI,:) = VdataToPlot{z}{ROI}{trialType}{trial};
                        VAVdataToPlot1{z}{trialType}{trial} = nanmean(VAVdataToPlot1_array{z}{trialType}{trial},1);

                        VAVdataToPlot2_array{trialType}{trial}(z,:) = VAVdataToPlot1{z}{trialType}{trial};
                        VAVdataToPlot{trialType}{trial} = nanmean(VAVdataToPlot2_array{trialType}{trial},1);
                     end 
                end 
            end 
        end 
    end 
end 

%average across trials
AVarray = cell(1,size(VAVdataToPlot,2));
AVdata = cell(1,size(VAVdataToPlot,2));
for trialType = 1:size(VAVdataToPlot,2)   
    if isempty(VAVdataToPlot{trialType}) == 0                  
        for trial = 1:length(VAVdataToPlot{trialType})  
            if isempty(VAVdataToPlot{trialType}{trial}) == 0                  
                AVarray{trialType}(trial,:) = VAVdataToPlot{trialType}{trial};
            end 
        end 
        SFWT6_ROI4av{trialType} = nanmean(AVarray{trialType},1);
     end 
end 

clear VAVdataToPlot1_array VAVdataToPlot1 VAVdataToPlot2_array VAVdataToPlot

%% average across imaging FOVs 

for trialType = 1:size(SFWT6_ROI3av,2)
    SFWT6av{trialType} = (SFWT6_ROI3av{trialType} + SFWT6_ROI4av{trialType})/2;
%     SF57av{trialType} = (SF57_ROI1av{trialType} + SF57_ROI2av{trialType})/2;
%     SF56av{trialType} = (SF56_ROI1av{trialType} + SF56_ROI2av{trialType})/2;
%     SF53av{trialType} = (SF53_ROI1av{trialType} + SF53_ROI2av{trialType})/2;
end 
% SF64av = SF64_ROI1av;

%% resample averaged individual mouse data 
for trialType = 1:size(SF63_ROI1av,2)
    if trialType == 1 || trialType == 3 
        goalLen = length(SF64av{1});
%         RSF58av{trialType} = resample(SF58av{trialType},goalLen,length(SF58av{trialType}));
        RSF63av{trialType} = resample(SF63av{trialType},goalLen,length(SF63av{trialType}));
%         RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
    elseif trialType == 2 || trialType == 4 
        goalLen = length(SF64av{4});
%         RSF58av{trialType} = resample(SF58av{trialType},goalLen,length(SF58av{trialType}));
        RSF63av{trialType} = resample(SF63av{trialType},goalLen,length(SF63av{trialType}));
%         RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
    end 
end 
RSF63av = SF63av;
RSF64av = SF64av;

 %% put all data into the same array 
% for trialType = 1:size(SF63av,2)
% %     if isempty(SF56av{trialType}) == 0 
%         miceData{trialType}(1,:) = RSF63av{trialType};
%         miceData{trialType}(2,:) = RSF64av{trialType};
% %         miceData{trialType}(3,:) = SF56av{trialType};
% %         miceData{trialType}(4,:) = RSF53av{trialType};
% %     end 
% end 
% %% get var and average across mice 
% for trialType = 1:size(SF61av,2)
%     varMiceData{trialType} = nanvar(miceData{trialType});
%     avMiceData{trialType} = nanmean(miceData{trialType});
% end 
miceData = SFWT6av;
avMiceData = SFWT6av;
%% smooth data if you want 
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ == 1 
    filtTime = input('How many seconds do you want to smooth your data by? '); 
    
    for trialType = 1:size(avMiceData,2)   
        if isempty(avMiceData{trialType}) == 0                  
            [VfiltD] = MovMeanSmoothData(miceData{trialType},filtTime,FPS);
%             [VfiltV] = MovMeanSmoothData(varMiceData{trialType},filtTime,FPS);
            VfiltData{trialType} = VfiltD;   
%             VfiltVar{trialType} = VfiltV; 
        end 
    end
     
elseif smoothQ == 0

    for trialType = 1:size(avMiceData,2)   
        if isempty(avMiceData{trialType}) == 0                           
%             VfiltData{trialType} = avMiceData{trialType};   
            VfiltData{trialType} = miceData{trialType};  
%             VfiltVar{trialType} = varMiceData{trialType};
        end 
    end
    
end 

%% plot 
dataMin = input("data Y axis MIN: ");
dataMax = input("data Y axis MAX: ");
FPSstack = FPS;
baselineEndFrame = round(20*(FPSstack));

for trialType = 4%1:size(miceData,2)  
    if isempty(miceData{trialType}) == 0

        figure;
        %set time in x axis            
        if trialType == 1 || trialType == 3 
            Frames = size(miceData{trialType}(1,:),2);                
            Frames_pre_stim_start = -((Frames-1)/2); 
            Frames_post_stim_start = (Frames-1)/2; 
            sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+2);
            FrameVals = round((1:FPSstack*2:Frames)-1); 
        elseif trialType == 2 || trialType == 4 
            Frames = size(miceData{trialType}(1,:),2);
            Frames_pre_stim_start = -((Frames-1)/2); 
            Frames_post_stim_start = (Frames-1)/2; 
            sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+10);
            FrameVals = round((1:FPSstack*2:Frames)-1); 
        end 
        ax=gca;
        plot(VfiltData{trialType},'k','LineWidth',2)
%         plot(VfiltData{trialType},'k')
        hold all;     

        %varargout = boundedline(1:size(VfiltData{trialType},2),VfiltData{trialType},VfiltVar{trialType},'k','transparency', 0.3,'alpha');                                                                             

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
            title(sprintf('Vessel data across mice smoothed by %d seconds.',filtTime));
        elseif smoothQ == 0
            title("Raw vessel data across mice.");
        end 

    end                       
end


