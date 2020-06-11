function [RdataToPlot,RdataToPlot2] = resampleVWdata(dataToPlot,dataToPlot2)

%% determine what dataset needs resampling 
% d1_len{1} = length(dataToPlot{1}{1}{1}{2});
% d1_len{2} = length(dataToPlot{1}{1}{2}{2});
% d1_len{3} = length(dataToPlot{1}{1}{3}{2});
d1_len{4} = length(dataToPlot{1}{1}{4}{2});

% d2_len{1} = length(dataToPlot2{1}{1}{1}{2});
% d2_len{2} = length(dataToPlot2{1}{1}{2}{2});
% d2_len{3} = length(dataToPlot2{1}{1}{3}{2});
d2_len{4} = length(dataToPlot2{1}{1}{4}{2});

if d1_len{4} > d2_len{4}
    dlenQ = 1;
elseif d1_len{4} < d2_len{4}
    dlenQ = 2;
end 

%% resample data 

if dlenQ == 1
    resamped_dataToPlot = cell(1,length(dataToPlot));
    for z = 1:length(dataToPlot)
        for VROI = 1:size(dataToPlot{z},2)
            for trialType = 1:size(dataToPlot{z}{VROI},2)
                if isempty(dataToPlot{z}{VROI}{trialType}) == 0 
                    for trial = 1:size(dataToPlot{z}{VROI}{trialType},2)
                        if isempty(dataToPlot{z}{VROI}{trialType}{trial}) == 0 
                            resamped_dataToPlot{z}{VROI}{trialType}{trial} = resample(dataToPlot{z}{VROI}{trialType}{trial},d2_len{trialType},d1_len{trialType});
                        end 
                    end 
                end 
            end 
        end 
    end 
    
   
    RdataToPlot = resamped_dataToPlot;
    RdataToPlot2 = dataToPlot2;
    
elseif dlenQ == 2
    resamped_dataToPlot2 = cell(1,length(dataToPlot2));
    for z = 1:length(dataToPlot2)
        for VROI = 1:size(dataToPlot2{z},2)
            for trialType = 1:size(dataToPlot2{z}{VROI},2)
                if isempty(dataToPlot2{z}{VROI}{trialType}) == 0 
                    for trial = 1:size(dataToPlot2{z}{VROI}{trialType},2)
                        if isempty(dataToPlot2{z}{VROI}{trialType}{trial}) == 0 
                            resamped_dataToPlot2{z}{VROI}{trialType}{trial} = resample(dataToPlot2{z}{VROI}{trialType}{trial},d1_len{trialType},d2_len{trialType});
                        end 
                    end 
                end 
            end 
        end 
    end 
    RdataToPlot = dataToPlot;
    RdataToPlot2 = resamped_dataToPlot2;
end 

end 