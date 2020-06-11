function plot_wheelAVdata(wheelDataToPlot,userInput,FPS,numZplanes,sec_before_stim_start)

%% smooth data if you want 
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ == 1 
    UIr = size(userInput,1)+1;
    filtTime = input('How many seconds do you want to smooth your data by? '); userInput(UIr,1) = ("How many seconds do you want to smooth your data by? "); userInput(UIr,2) = (filtTime); UIr = UIr+1;
    
    WfiltData = cell(1,size(wheelDataToPlot,2));
    for trialType = 1:size(wheelDataToPlot,2)             
        if isempty(wheelDataToPlot{trialType}) == 0 
            for trial = 1:size(wheelDataToPlot{trialType},2) 
                if isempty(wheelDataToPlot{trialType}{trial})==0
                    [WfiltD] = MovMeanSmoothData(wheelDataToPlot{trialType}{trial},filtTime,FPS);             
                    WfiltData{trialType}{trial} = WfiltD;
                end 
            end 
        end 
    end 
    
elseif smoothQ == 0 
   
    WfiltData = cell(1,size(wheelDataToPlot,2));
    for trialType = 1:size(wheelDataToPlot,2)             
        if isempty(wheelDataToPlot{trialType}) == 0 
            for trial = 1:size(wheelDataToPlot{trialType},2)                         
                WfiltData{trialType}{trial} = wheelDataToPlot{trialType}{trial};
            end 
        end 
    end 
    
end 

% average across trials
WAVarray = cell(1,size(wheelDataToPlot,2));
WAVdata = cell(1,size(wheelDataToPlot,2));
WSEMdata = cell(1,size(wheelDataToPlot,2));
for trialType = 1:size(wheelDataToPlot,2)
    if isempty(wheelDataToPlot{trialType}) == 0 
        for trial = 1:size(wheelDataToPlot{trialType},2)
            if isempty(WfiltData{trialType}{trial})==0
                WAVarray{trialType}(trial,:) = WfiltData{trialType}{trial};
            end 
        end 
        WAVdata{trialType}= nanmean(WAVarray{trialType},1);
        WSEMdata{trialType} = (nanstd(WAVarray{trialType},1))/(sqrt(size(WfiltData{trialType},2)));
    end      
end 


%% plot 
dataMin = input("data Y axis MIN: ");
dataMax = input("data Y axis MAX: ");
FPSstack = FPS/numZplanes;
baselineEndFrame = round(sec_before_stim_start*(FPSstack));

%for count = length(CfiltData)
    for trialType = 1:size(WAVdata,2)  
        if isempty(WAVdata{trialType}) == 0

            figure;
            %set time in x axis            
            if trialType == 1 || trialType == 3 
                Frames = size(wheelDataToPlot{trialType}{2},2);                
                Frames_pre_stim_start = -((Frames-1)/2); 
                Frames_post_stim_start = (Frames-1)/2; 
                sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+2);
                FrameVals = round((1:FPSstack*2:Frames)-1); 
            elseif trialType == 2 || trialType == 4 
                Frames = size(wheelDataToPlot{trialType}{2},2);
                Frames_pre_stim_start = -((Frames-1)/2); 
                Frames_post_stim_start = (Frames-1)/2; 
                sec_TimeVals = floor(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+11);
                FrameVals = round((1:FPSstack*2:Frames)-1); 
            end 
            ax=gca;
            plot(WAVdata{trialType},'k')
            hold all;     

            varargout = boundedline(1:size(WAVdata{trialType},2),WAVdata{trialType},WSEMdata{trialType},'k','transparency', 0.3,'alpha');                                                                             

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

            ylim([dataMin dataMax]);

            if smoothQ == 1 
                title(sprintf('Wheel data smoothed by %d seconds.',filtTime));
            elseif smoothQ == 0
                title("Raw Wheel Data.");
            end 

        end                       
    end




end 