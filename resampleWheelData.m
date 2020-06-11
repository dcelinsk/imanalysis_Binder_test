function [RWdataToPlot,RWdataToPlot2] = resampleWheelData(WdataToPlot,WdataToPlot2)

%% determine what dataset needs resampling 
% d1_len{1} = length(WdataToPlot{1}{2});
% d1_len{2} = length(WdataToPlot{2}{2});
% d1_len{3} = length(WdataToPlot{3}{2});
d1_len{4} = length(WdataToPlot{4}{2});
% 
% d2_len{1} = length(WdataToPlot2{1}{2});
% d2_len{2} = length(WdataToPlot2{2}{2});
% d2_len{3} = length(WdataToPlot2{3}{2});
d2_len{4} = length(WdataToPlot2{4}{2});

if d1_len{4} > d2_len{4}
    dlenQ = 1;
elseif d1_len{4} < d2_len{4}
    dlenQ = 2;
end 

%% resample data 

if dlenQ == 1
    for trialType = 1:size(WdataToPlot,2)
        if isempty(WdataToPlot{trialType}) == 0                
            for trial = 1:size(WdataToPlot{trialType},2)                  
                if isempty(WdataToPlot{trialType}{trial}) == 0                   
                    resamped_WdataToPlot{trialType}{trial} = resample(WdataToPlot{trialType}{trial},d2_len{trialType},d1_len{trialType});                  
                end 
            end 
        end 
    end     

    RWdataToPlot = resamped_WdataToPlot;
    RWdataToPlot2 = WdataToPlot2;
    
elseif dlenQ == 2
    resamped_WdataToPlot2 = cell(1,length(WdataToPlot2));
    for trialType = 1:size(WdataToPlot2,2)
        if isempty(WdataToPlot2{trialType}) == 0                
            for trial = 1:size(WdataToPlot2{trialType},2)                  
                if isempty(WdataToPlot2{trialType}{trial}) == 0 
                    resamped_WdataToPlot2{trialType}{trial} = resample(WdataToPlot2{trialType}{trial},d1_len{trialType},d2_len{trialType});
                end 
            end 
        end 
    end     
    
    RWdataToPlot = WdataToPlot;
    RWdataToPlot2 = resamped_WdataToPlot2;
end 

end 