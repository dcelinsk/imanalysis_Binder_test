%% get just the BBB data
temp = matfile('SF64_20200222_70RhoB_ROI1_DAca.mat');
SF64_ROI1_Cdata = temp.dataToPlot;

%% resample 
% [RSF58_ROI1_Cdata,RSF58_ROI2_Cdata] = resampleCadata(SF58_ROI1_Cdata,SF58_ROI2_Cdata);
[RSF63_ROI1_Cdata,RSF63_ROI2_Cdata] = resampleCadata(SF63_ROI1_Cdata,SF63_ROI2_Cdata);
% [RSF56_ROI1_Bdata,RSF56_ROI2_Bdata] = resampleBBBdata(SF56_ROI1_Bdata,SF56_ROI2_Bdata);
% [RSF53_ROI1_Bdata,RSF53_ROI2_Bdata] = resampleBBBdata(SF53_ROI1_Bdata,SF53_ROI2_Bdata);

RSF64_ROI1_Cdata = SF64_ROI1_Cdata;


%% average across planes in Z, ROIs, and trials
CdataToPlot = RSF64_ROI1_Cdata;

count = 1;
CAVdataToPlot1_array = cell(1,length(CdataToPlot));
CAVdataToPlot1 = cell(1,length(CdataToPlot));
CAVdataToPlot2_array = cell(1,length(CdataToPlot));
for posCell = 1:length(CdataToPlot)
    if isempty(CdataToPlot{posCell}) == 0             
        for trialType = 1:size(CdataToPlot{posCell},2)
            for z = 1:size(CdataToPlot{posCell},1)           
                for trial = 1:size(CdataToPlot{posCell}{z,trialType},2)  
                    if isempty(CdataToPlot{posCell}{z,trialType}{trial}) == 0
                    
                        CAVdataToPlot1_array{posCell}{trialType}{trial}(z,:) = CdataToPlot{posCell}{z,trialType}{trial};
                        CAVdataToPlot1{posCell}{trialType}{trial} = nanmean(CAVdataToPlot1_array{posCell}{trialType}{trial},1);


                        CAVdataToPlot2_array{posCell}{trialType}(trial,:) = CAVdataToPlot1{posCell}{trialType}{trial};                    
                        CAVdataToPlot{count}{trialType} = nanmean(CAVdataToPlot2_array{posCell}{trialType},1);
                    end 
                end 
            end 
        end        
        count = count + 1;
    end 
end 

%average across cells 
AVarray = cell(1,length(CAVdataToPlot));
AVdata = cell(1,length(CAVdataToPlot));
for count = 1:length(CAVdataToPlot)
    for trialType = 1:size(CAVdataToPlot{1},2)
        if isempty(CAVdataToPlot{count}{trialType}) == 0 
            AVarray{trialType}(count,:) = CAVdataToPlot{count}{trialType};
            AVdata{trialType} = nanmean(AVarray{trialType},1);
        end 
    end     
end 

SF64_ROI1av = AVdata;

clear CAVdataToPlot1_array CAVdataToPlot1 CAVdataToPlot2_array CAVdataToPlot AVarray AVdata

%% average across imaging FOVs 

for trialType = 1:4%size(BdataToPlot{2},2)
%     SF58av{trialType} = SF58_ROI1av{trialType};%(SF58_ROI1av{trialType} + SF58_ROI2av{trialType})/2;
    SF63av{trialType} = (SF63_ROI1av{trialType} + SF63_ROI2av{trialType})/2;
%     SF56av{trialType} = SF56_ROI2av{trialType};%(SF56_ROI1av{trialType} + SF56_ROI2av{trialType})/2;
%     SF53av{trialType} = SF53_ROI1av{trialType};%(SF53_ROI1av{trialType} + SF53_ROI2av{trialType})/2;
end 
SF64av = SF64_ROI1av;
% SF58av = SF58_ROI1av;

%% resample averaged individual mouse data 
% for trialType = 1:4
%     if isempty(SF57av{trialType}) == 0 
%         if trialType == 1 || trialType == 3 
%             goalLen = length(SF58av{1});
%             RSF56av{trialType} = resample(SF56av{trialType},goalLen,length(SF56av{trialType}));
%             RSF57av{trialType} = resample(SF57av{trialType},goalLen,length(SF57av{trialType}));
% %             RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
%         elseif trialType == 2 || trialType == 4 
%             goalLen = length(SF58av{2});
%             RSF56av{trialType} = resample(SF56av{trialType},goalLen,length(SF56av{trialType}));
%             RSF57av{trialType} = resample(SF57av{trialType},goalLen,length(SF57av{trialType}));
% %             RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
%         end 
%     end 
% end 
RSF63av = SF63av;
RSF64av = SF64av;
%% put all data into the same array 
for trialType = 1:4
    if isempty(SF64av{trialType}) == 0 
%         miceData{trialType}(1,:) = RSF58av{trialType};
        miceData{trialType}(1,:) = RSF63av{trialType};
        miceData{trialType}(2,:) = RSF64av{trialType};
%         miceData{trialType}(4,:) = RSF53av{trialType};
    end 
end 

%% get SEM and average across mice 
for trialType = 1:size(SF64_ROI1av,2)
    varMiceData{trialType} = nanvar(miceData{trialType});
    avMiceData{trialType} = nanmean(miceData{trialType});
end 

%% smooth data if you want 
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ == 1 
    filtTime = input('How many seconds do you want to smooth your data by? '); 
    
    for trialType = 1:size(avMiceData,2)   
        if isempty(avMiceData{trialType}) == 0                  
            [VfiltD] = MovMeanSmoothData(avMiceData{trialType},filtTime,FPS);
            [VfiltV] = MovMeanSmoothData(varMiceData{trialType},filtTime,FPS);
            VfiltData{trialType} = VfiltD;   
            VfiltVar{trialType} = VfiltV; 
        end 
    end
     
elseif smoothQ == 0

    for trialType = 1:size(avMiceData,2)   
        if isempty(avMiceData{trialType}) == 0                           
            VfiltData{trialType} = avMiceData{trialType};   
            VfiltVar{trialType} = varMiceData{trialType};
        end 
    end
    
end 

%% plot 
dataMin = input("data Y axis MIN: ");
dataMax = input("data Y axis MAX: ");
FPSstack = FPS/3;
%FPSstack = FPS;
baselineEndFrame = round(20*(FPSstack));



for trialType = 1:size(avMiceData,2)  
    if isempty(avMiceData{trialType}) == 0

        figure;
        %set time in x axis            
        if trialType == 1 || trialType == 3 
            Frames = size(avMiceData{trialType},2);                
            Frames_pre_stim_start = -((Frames-1)/2); 
            Frames_post_stim_start = (Frames-1)/2; 
            sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+2);
            FrameVals = round((1:FPSstack*2:Frames)-1); 
        elseif trialType == 2 || trialType == 4 
            Frames = size(avMiceData{trialType},2);
            Frames_pre_stim_start = -((Frames-1)/2); 
            Frames_post_stim_start = (Frames-1)/2; 
            sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+11);
            FrameVals = round((1:FPSstack*2:Frames)-1); 
        end 
        ax=gca;
        plot(VfiltData{trialType},'r')
        hold all;     

        varargout = boundedline(1:size(VfiltData{trialType},2),VfiltData{trialType},VfiltVar{trialType},'b','transparency', 0.3,'alpha');                                                                             

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
            title(sprintf('DA terminal calcium data across mice smoothed by %d seconds.',filtTime));
        elseif smoothQ == 0
            title("Raw DA terminal calcium data across mice.");
        end 

    end                       
end


