%% get just the data you need 
temp = matfile('70RhoB_DAT-GCaMP-SF63-64_BBB.mat');
BBBdata = temp.miceData;


%% put data in same cell array for simplicity moving forward 
data(1,1:4) = DAT_Chrimson_GCaMP_data;
data(2,1:4) = DAT_Chrimson_data;
data(1,1:4) = DAT_GCaMP_data;
data(4,1:4) = WT_data;

data(1,1:4) = DAT_Chrimson_GCaMP_70FITC_BBB_data;
data(2,1:4) = DAT_Chrimson_GCaMP_10FITC_BBB_data;
data(3,1:4) = DAT_Chrimson_GCaMP_3FITC_BBB_data;

%this separates chrimson from non chrimson mice 
data2(1,1:4) = R_data(1,1:4);
for trialType = 1:4
    data2{1,trialType}(5,:) = R_data{2,trialType}(1,:);
end 
data2(2,1:4) = R_data(3,1:4);
for trialType = 4
    data2{2,trialType}(3,:) = R_data{4,trialType}(1,:);
end 

%% resample miceData 

for dataType = 1:4
    for trialType = 1:4 
        if trialType == 1 || trialType == 3 
            goalLength = length(DAT_GCaMP_data{1});
        elseif trialType == 2 || trialType == 4 
            goalLength = length(DAT_GCaMP_data{4});  
        end 
        for mouse = 1:size(data{dataType,trialType},1)
            R_data{dataType,trialType}(mouse,:) = resample(data{dataType,trialType}(mouse,:),goalLength,length(data{dataType,trialType}));
        end 
    end 
end 

%% 
%  for dataType = 1:4
%     for trialType = 1:4
%         avData_0{dataType,trialType} = nanmean(R_data{dataType,trialType},1);
%         varData_0{dataType,trialType} = (nanstd(R_data{dataType,trialType},1))/(sqrt(size(S_data{dataType,trialType},1)));
%     end
%  end 

%% smooth data if you want 
% normalize data to baseline period - plot % change 
baselineEnd = ceil((FPS)*18);
for dataType =  1%:size(data,1)
    for trialType = 1:size(data,2)
        for mouse = 1:size(R_data{dataType,trialType},1)
%             N_data{dataType,trialType}(mouse,:) = (R_data{dataType,trialType}(mouse,:))/abs((mean(R_data{dataType,trialType}(mouse,1:baselineEnd))));
            N_data{dataType,trialType}(mouse,:) = (R_data{dataType,trialType}(mouse,:))/abs((mean(R_data{dataType,trialType}(mouse,1:baselineEnd))));

        end 
    end 
end 


%{
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ == 1 
    sigmaf = 30;
    for dataType = 1:4
        for trialType = 1:4
            for mouse = 1:size(N_data{dataType,trialType},1)
                [sData] = FftRft(N_data{dataType,trialType}(mouse,:),FPS,sigmaf);
                S_data{dataType,trialType}(mouse,:) = sData;
            end 
        end
    end 
     
elseif smoothQ == 0
    S_data = N_data;
end
%}
smoothQ = input('Do you want to smooth your data? Yes = 1. No = 0. ');

if smoothQ ==  1
    filtTime = input('How many seconds do you want to smooth your data by? '); 
    for dataType = 1:size(data,1)
        for trialType = 1:size(data,2)
            for mouse = 1:size(R_data{dataType,trialType},1)
                [sData] = MovMeanSmoothData(R_data{dataType,trialType}(mouse,:),filtTime,FPS);
                S_data{dataType,trialType}(mouse,:) = sData;
            end 
%             [sData] = MovMeanSmoothData(avData_0{dataType,trialType}(mouse,:),filtTime,FPS);
%             S_data{dataType,trialType} = sData;
        end
    end 
     
elseif smoothQ == 0
    S_data = N_data;
end

% 
 % get var and average across mice 
 for dataType = 1:size(data,1)
    for trialType = 1:size(data,2)
        avData{dataType,trialType} = nanmean(S_data{dataType,trialType},1);
        varData{dataType,trialType} = (nanstd(S_data{dataType,trialType},1))/(sqrt(size(S_data{dataType,trialType},1)));
    end
 end 
 
% avData = S_data;

%%
dataMin = input("data Y axis MIN: ");
dataMax = input("data Y axis MAX: ");
  
%% plot

FPSstack = FPS;
baselineEndFrame = ceil((FPS)*20);

for trialType = 2
    figure;
    hold all;

    %set time in x axis            
    if trialType == 1 || trialType == 3 
        Frames = size(avData{1,trialType}(1,:),2);                
        Frames_pre_stim_start = -((Frames-1)/2); 
        Frames_post_stim_start = (Frames-1)/2; 
        sec_TimeVals = ceil(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+0);
        FrameVals = round((1:FPSstack*2:Frames)-1); 
    elseif trialType == 2 || trialType == 4 
        Frames = size(avData{1,trialType}(1,:),2);
        Frames_pre_stim_start = -((Frames-1)/2); 
        Frames_post_stim_start = (Frames-1)/2; 
        sec_TimeVals = ceil(((Frames_pre_stim_start:FPSstack*2:Frames_post_stim_start)/FPSstack)+9);
        FrameVals = round((1:FPSstack*2:Frames)-1); 
    end 


    ax=gca;
    plot(avData{1,trialType},'k','LineWidth',2)  
    plot(BavData{1,trialType},'r','LineWidth',2)   
    
%     plot(VavData{2,trialType},'r','LineWidth',2)  
%     plot(avData{4,trialType},'k','LineWidth',2)  

%       plot(RdiffData,'k','LineWidth',2)
%       plot(avData{1,trialType},'r','LineWidth',2)
      
% figure ;
%     varargout = boundedline(1:size(RCAdata,2),RCAdata,RCAvar,'g','transparency', 0.3,'alpha');                                                                             
%     varargout = boundedline(1:size(avData{1,trialType},2),avData{1,trialType},varData{1,trialType},'k','transparency', 0.3,'alpha'); 

    varargout = boundedline(1:size(BavData{1,trialType},2),BavData{1,trialType},BvarData{1,trialType},'r','transparency', 0.3,'alpha');                                                                             
    varargout = boundedline(1:size(avData{1,trialType},2),avData{1,trialType},varData{1,trialType},'k','transparency', 0.3,'alpha'); 
%     varargout = boundedline(1:size(VavData{2,trialType},2),VavData{2,trialType},VvarData{2,trialType},'r','transparency', 0.3,'alpha');                                                                             


    ax.XTick = FrameVals;
    ax.XTickLabel = sec_TimeVals;
    ax.FontSize = 20;
    if trialType == 1 
        plot([(baselineEndFrame+((FPSstack)*4))-0.3 (baselineEndFrame+((FPSstack)*4))-0.3], [-10000000 10000000], 'Color',[0.3,0.7,0.9],'LineWidth',3)
        plot([(baselineEndFrame+((FPSstack)*2))-0 (baselineEndFrame+((FPSstack)*2))-0], [-10000000 10000000], 'b','LineWidth',3)
        %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*2)) round(baselineEndFrame+((FPSstack/numZplanes)*2)) baselineEndFrame],[-5000 -5000 5000 5000],'b')
        %alpha(0.4)   
        plot([baselineEndFrame baselineEndFrame], [-10000000 10000000], 'b','LineWidth',3) 
        plot([baselineEnd baselineEnd], [-10000000 10000000], 'Color',[0.3,0.7,0.9],'LineWidth',3) 
    elseif trialType == 3 
%         plot([(baselineEndFrame+((FPSstack)*4))-0.3 (baselineEndFrame+((FPSstack)*4))-0.3 ], [-10000000 10000000], 'Color',[0.8,0.5,0.5],'LineWidth',3)
        plot([round(baselineEndFrame+((FPSstack)*2))-1 round(baselineEndFrame+((FPSstack)*2))-1], [-10000000 10000000], 'r','LineWidth',3)
        %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*2)) round(baselineEndFrame+((FPSstack/numZplanes)*2)) baselineEndFrame],[-5000 -5000 5000 5000],'r')
        %alpha(0.4)     
        plot([baselineEndFrame baselineEndFrame], [-10000000 100000000], 'r','LineWidth',3) 
%         plot([baselineEnd baselineEnd], [-10000000 10000000], 'Color',[0.8,0.5,0.5],'LineWidth',3) 
    elseif trialType == 2 
%         plot([(baselineEndFrame+((FPSstack)*22))-0  (baselineEndFrame+((FPSstack)*22))-0], [-10000000 10000000], 'Color',[0.3,0.7,0.9],'LineWidth',3)
        plot([round(baselineEndFrame+((FPSstack)*20))-0 round(baselineEndFrame+((FPSstack)*20))-0], [-100000000 1000000], 'b','LineWidth',3)
        %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*20)) round(baselineEndFrame+((FPSstack/numZplanes)*20)) baselineEndFrame],[-5000 -5000 5000 5000],'b')
        %alpha(0.4)   
        plot([baselineEndFrame baselineEndFrame], [-10000000 10000000], 'b','LineWidth',3) 
%         plot([baselineEnd baselineEnd], [-10000000 10000000], 'Color',[0.3,0.7,0.9],'LineWidth',3) 
    elseif trialType == 4 
%         plot([(baselineEndFrame+((FPSstack)*22))-0 (baselineEndFrame+((FPSstack)*22))-0], [-10000000 10000000], 'Color',[0.8,0.5,0.5],'LineWidth',3) 
        plot([round(baselineEndFrame+((FPSstack)*20))-0.5 round(baselineEndFrame+((FPSstack)*20))-0.5], [-10000000 10000000], 'r','LineWidth',3)
        %patch([baselineEndFrame round(baselineEndFrame+((FPSstack/numZplanes)*20)) round(baselineEndFrame+((FPSstack/numZplanes)*20)) baselineEndFrame],[-5000 -5000 5000 5000],'r')
        %alpha(0.4)  
        plot([baselineEndFrame baselineEndFrame], [-10000000 10000000], 'r','LineWidth',3) 
%         plot([baselineEnd-0.4 baselineEnd-0.4], [-10000000 10000000], 'Color',[0.8,0.5,0.5],'LineWidth',3) 
    end
   
%     legend('DAT-Chrimson-GCaMP (N = 4)','DAT-Chrimson (N = 1)','DAT-GCaMP (N = 2)')
%     legend('DAT-Chrimson-GCaMP (N = 4)','DAT-Chrimson (N = 1)','WT (N = 1)')
%     legend('DAT-GCaMP (N = 2)')
%     legend('DAT-Chrimson-GCaMP (N = 4)','DAT-Chrimson (N = 1)','DAT-GCaMP (N = 2)','WT (N = 1)')
%     legend('DAT-GCaMP (N = 2)')
%     legend('DAT-Chrimson-GCaMP (N = 4)','DAT-GCaMP (N = 2)')
%     legend('DAT-Chrimson-GCaMP (N = 4)')
%     legend('DAT-Chrimson (N = 1)','DAT-GCaMP (N = 2)')
%     legend('DAT-Chrimson (N = 1)','DAT-GCaMP (N = 2)','WT (N = 1)')
%     legend('DAT-Chrimson (N = 1)','WT (N = 1)')
%     legend('70 kD FITC (N = 4)','10 kD FITC (N = 4)')%,'3 kD FITC (N = 3)')
%     legend('DAT-Chrimson-GCaMP (N = 4)','No Chrimson (N = 3)')
%     legend('No Chrimson (N = 3)')
%     legend('BBB data','Vessel width')
%     legend('DAT-Chrimson-GCaMP (N = 4)','DAT-Chrimson (N = 1)')
%     legend('DAT-Chrimson-GCaMP (N = 5)')%,'DAT-GCaMP (N = 2)')
%     legend('GCaMP data','DAT-Chrimson-GCaMP - No-Chrimson Data','Vessel Width')
%     legend('No Chrimson (N = 3)')
    legend('BBB data','Vessel width')

    ylim([dataMin dataMax]);
    if trialType == 1 || trialType == 3 
        xlim([1 length(avData{1,trialType})])
    elseif trialType == 2 || trialType == 4 
        xlim([1 length(avData{1,trialType})])
    end 
    xlabel('time (s)');
%     ylabel('percent change from baseline')
    ylabel('z-scored data')

    if smoothQ == 1 
        title(sprintf('No-Chrimson (N = 2) data smoothed by %d seconds',filtTime));
    elseif smoothQ == 0
        title("Raw Vessel width across mice.");
    end 

                     
 end 