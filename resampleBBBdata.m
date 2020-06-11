function [RBdataToPlot,RBdataToPlot2] = resampleBBBdata(BdataToPlot,BdataToPlot2)

%% determine what dataset needs resampling 
% d1_len{1} = length(BdataToPlot{1}{1}{1}{1});
% d1_len{2} = length(BdataToPlot{1}{2}{1}{1});
% d1_len{3} = length(BdataToPlot{1}{3}{1}{1});
d1_len{4} = length(BdataToPlot{1}{4}{1}{1});

% d2_len{1} = length(BdataToPlot2{1}{1}{1}{1});
% d2_len{2} = length(BdataToPlot2{1}{2}{1}{1});
% d2_len{3} = length(BdataToPlot2{1}{3}{1}{1});
d2_len{4} = length(BdataToPlot2{1}{4}{1}{1});

if d1_len{4} > d2_len{4}
    dlenQ = 1;
elseif d1_len{4} < d2_len{4}
    dlenQ = 2;
end 

%% resample data 

if dlenQ == 1
    resamped_BdataToPlot = cell(1,length(BdataToPlot));
    for z = 1:length(BdataToPlot)
        for trialType = 1:size(BdataToPlot{z},2)
            if isempty(BdataToPlot{z}{trialType}) == 0                
                for trial = 1:size(BdataToPlot{z}{trialType},2)                  
                    if isempty(BdataToPlot{z}{trialType}{trial}) == 0 
                        for ROI = 1:size(BdataToPlot{z}{trialType}{trial},2) 
                            resamped_BdataToPlot{z}{trialType}{trial}{ROI} = resample(BdataToPlot{z}{trialType}{trial}{ROI},d2_len{trialType},d1_len{trialType});
                        end 
                    end 
                end 
            end 
        end     
    end 
    
   
    RBdataToPlot = resamped_BdataToPlot;
    RBdataToPlot2 = BdataToPlot2;
    
elseif dlenQ == 2
    resamped_BdataToPlot2 = cell(1,length(BdataToPlot2));
    for z = 1:length(BdataToPlot2)
        for trialType = 1:size(BdataToPlot2{z},2)
            if isempty(BdataToPlot2{z}{trialType}) == 0                
                for trial = 1:size(BdataToPlot2{z}{trialType},2)                  
                    if isempty(BdataToPlot2{z}{trialType}{trial}) == 0 
                        for ROI = 1:size(BdataToPlot2{z}{trialType}{trial},2) 
                            resamped_BdataToPlot2{z}{trialType}{trial}{ROI} = resample(BdataToPlot2{z}{trialType}{trial}{ROI},d1_len{trialType},d2_len{trialType});
                        end 
                    end 
                end 
            end 
        end     
    end 
    RBdataToPlot = BdataToPlot;
    RBdataToPlot2 = resamped_BdataToPlot2;
end 

end 