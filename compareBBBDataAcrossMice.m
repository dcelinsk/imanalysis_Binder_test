%% get just the BBB data
temp = matfile('SFWT6_20190607_ROI4_BBB.mat');
SFWT6_ROI4_Bdata = temp.BdataToPlot;

%% resample 
[RSFWT6_ROI3_Bdata,RSFWT6_ROI4_Bdata] = resampleBBBdata(SFWT6_ROI3_Bdata,SFWT6_ROI4_Bdata);
% [RSF57_ROI1_Bdata,RSF57_ROI2_Bdata] = resampleBBBdata(SF57_ROI1_Bdata,SF57_ROI2_Bdata);
% [RSF56_ROI1_Bdata,RSF56_ROI2_Bdata] = resampleBBBdata(SF56_ROI1_Bdata,SF56_ROI2_Bdata);
% [RSF53_ROI1_Bdata,RSF53_ROI2_Bdata] = resampleBBBdata(SF53_ROI1_Bdata,SF53_ROI2_Bdata);

%% average across planes in Z, ROIs, and trials
% BdataToPlot = RSFWT6_ROI4_Bdata;
BdataToPlot = SF_57_70FITC_ROI2_Bdata;

for Z = 1:length(BdataToPlot)
    for trialType = 1:size(BdataToPlot{Z},2)        
        if isempty(BdataToPlot{Z}{trialType}) == 0 
            for trial = 1:size(BdataToPlot{Z}{trialType},2)
                for VROI = 1:size(BdataToPlot{Z}{trialType}{trial},2)
                    BAVdataToPlot1_array{Z}{trialType}{trial}(VROI,:) = BdataToPlot{Z}{trialType}{trial}{VROI};
                    BAVdataToPlot1{Z}{trialType}{trial} = nanmean(BAVdataToPlot1_array{Z}{trialType}{trial},1);

                    BAVdataToPlot2_array{trialType}{trial}(Z,:) = BAVdataToPlot1{Z}{trialType}{trial};
                    BAVdataToPlot{trialType}{trial} = nanmean(BAVdataToPlot2_array{trialType}{trial},1);
                end 
            end 
        end 
        
    end 
end 

%average across trials
AVarray = cell(1,length(BAVdataToPlot));
AVdata = cell(1,length(BAVdataToPlot));
for trialType = 1:size(BAVdataToPlot,2)
    if isempty(BAVdataToPlot{trialType}) == 0 
        for trial = 1:size(BAVdataToPlot{trialType},2)
            AVarray{trialType}(trial,:) = BAVdataToPlot{trialType}{trial};
            AVdata{trialType} = nanmean(AVarray{trialType},1);
        end 
    end   
end    

SF57_ROI2av = AVdata;

% clear BAVdataToPlot1_array BAVdataToPlot1 BAVdataToPlot2_array BAVdataToPlot

%% average across imaging FOVs 

for trialType = 4%1:size(SFWT6_ROI1av,2)
    SFWT6av{trialType} = (SFWT6_ROI3av{trialType} + SFWT6_ROI4av{trialType})/2;
%     SF57av{trialType} = (SF57_ROI1av{trialType} + SF57_ROI2av{trialType})/2;
%     SF56av{trialType} = SF56_ROI2av{trialType};%(SF56_ROI1av{trialType} + SF56_ROI2av{trialType})/2;
%     SF53av{trialType} = SF53_ROI1av{trialType};%(SF53_ROI1av{trialType} + SF53_ROI2av{trialType})/2;
end 
SF56av = SF56_ROI2av;
SF58av = SF58_ROI2av;

%% resample averaged individual mouse data 
for trialType = 1:size(SF58_ROI2av,2)
    if trialType == 1 || trialType == 3 
        goalLen = length(SF58av{1});
%         RSF58av{trialType} = resample(SF58av{trialType},goalLen,length(SF58av{trialType}));
        RSF56av{trialType} = resample(SF56av{trialType},goalLen,length(SF56av{trialType}));
%         RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
    elseif trialType == 2 || trialType == 4 
        goalLen = length(SF58av{4});
%         RSF58av{trialType} = resample(SF58av{trialType},goalLen,length(SF58av{trialType}));
        RSF56av{trialType} = resample(SF56av{trialType},goalLen,length(SF56av{trialType}));
%         RSF53av{trialType} = resample(SF53av{trialType},goalLen,length(SF53av{trialType}));
    end 
end 
RSF58av = SF58av;
% RSFWT6av = SFWT6av;

%% put all data into the same array 
for trialType = 1:size(SF58_ROI2av,2)
    if isempty(RSF56av{trialType}) == 0 
        miceData{trialType}(1,:) = RSF56av{trialType};
        miceData{trialType}(2,:) = RSF58av{trialType};
%         miceData{trialType}(3,:) = RSF58av{trialType};
%         miceData{trialType}(4,:) = RSF53av{trialType};
    end 
end 


%% smooth data if you want 
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ == 1 
    filtTime = input('How many seconds do you want to smooth your data by? '); 
    
    for trialType = 1:size(miceData,2)   
        if isempty(miceData{trialType}) == 0                  
            [VfiltD] = MovMeanSmoothData(miceData{trialType},filtTime,FPS);
%             [BVfiltD] = MovMeanSmoothData(BmiceData{trialType},filtTime,FPS);
            [VfiltV] = MovMeanSmoothData(varMiceData{trialType},filtTime,FPS);
            VfiltData{trialType} = VfiltD;   
%             BVfiltData{trialType} = BVfiltD;  
            VfiltVar{trialType} = VfiltV; 
        end 
    end
     
elseif smoothQ == 0

    for trialType = 1:size(miceData,2)   
        if isempty(miceData{trialType}) == 0                           
            VfiltData{trialType} = miceData{trialType}; 
%             BVfiltData{trialType} = BmiceData{trialType}; 
            VfiltVar{trialType} = varMiceData{trialType};
        end 
    end
    
end 

% get SEM and average across mice 
for trialType = 1:size(VfiltData,2)
    
    avMiceData{trialType} = nanmean(VfiltData{trialType},1);
    
    varMiceData{trialType} = (nanstd(VfiltData{trialType},1)).^2;
    semMiceData{trialType} = ((nanstd(VfiltData{trialType},1)))/(sqrt(size(VfiltData{trialType},1)));
%     
%     BavMiceData{trialType} = nanmean(BVfiltData{trialType},1);
%     
%     BsemMiceData{trialType} = ((nanstd(BVfiltData{trialType},1)))/(sqrt(size(BVfiltData{trialType},1)));
end 

%  VSEMdata{trialType} = (nanstd(VAVarray{trialType},1))/(sqrt(size(VAVdataToPlot{trialType},2)));


%% plot 
% dataMin = input("data Y axis MIN: ");
% dataMax = input("data Y axis MAX: ");
FPSstack = FPS;%/3;
baselineEndFrame = round(20*(FPSstack));



for trialType = 1:size(miceData,2)  
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
        plot(avMiceData{trialType},'b','LineWidth',2)
        hold all;     
%         plot(BavMiceData{trialType},'r','LineWidth',2)

        varargout = boundedline(1:size(avMiceData{trialType},2),avMiceData{trialType},semMiceData{trialType},'b','transparency', 0.3,'alpha');                                                                             
%         varargout = boundedline(1:size(BavMiceData{trialType},2),BavMiceData{trialType},BsemMiceData{trialType},'r','transparency', 0.3,'alpha');                                                                             


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
%             legend('DA calcium','BBB data')
        ylim([dataMin dataMax]);

        if smoothQ == 1 
            title(sprintf('Data across mice smoothed by %d seconds.',filtTime));
        elseif smoothQ == 0
            title("Raw data across mice.");
        end 

    end                       
end


