function [RBdataToPlot,RBdataToPlot2] = resampleCadata(BdataToPlot,BdataToPlot2)

%% determine what dataset needs resampling 
for posCell = 1:length(BdataToPlot)
    if isempty(BdataToPlot{posCell}) == 0
        d1_len{1} = length(BdataToPlot{posCell}{1,1}{1});
        d1_len{2} = length(BdataToPlot{posCell}{1,2}{1});
        d1_len{3} = length(BdataToPlot{posCell}{1,3}{1});
        d1_len{4} = length(BdataToPlot{posCell}{1,4}{1});
    end 
end 

for posCell = 1:length(BdataToPlot2)
    if isempty(BdataToPlot2{posCell}) == 0
        d2_len{1} = length(BdataToPlot2{posCell}{1,1}{1});
        d2_len{2} = length(BdataToPlot2{posCell}{1,2}{1});
        d2_len{3} = length(BdataToPlot2{posCell}{1,3}{1});
%         d2_len{4} = length(BdataToPlot2{posCell}{1,4}{1});
    end 
end 

if d1_len{1} > d2_len{1}
    dlenQ = 1;
elseif d1_len{1} < d2_len{1}
    dlenQ = 2;
end 

%% resample data 

if dlenQ == 1
    resamped_BdataToPlot = cell(1,length(BdataToPlot));
    for posCell = 1:length(BdataToPlot)
        if isempty(BdataToPlot{posCell}) == 0             
            for trialType = 1:size(BdataToPlot{posCell},2)    
                if trialType == 1 || trialType == 3 
                    for z = 1:size(BdataToPlot{posCell},1)    
                         if isempty(BdataToPlot{posCell}{z,trialType}) == 0
                            for trial = 1:size(BdataToPlot{posCell}{z,trialType},2)      
                                if isempty(BdataToPlot{posCell}{z,trialType}{trial}) == 0 
                                    resamped_BdataToPlot{posCell}{z,trialType}{trial} = resample(BdataToPlot{posCell}{z,trialType}{trial},d2_len{1},d1_len{trialType});                                                                    
                                end 
                            end 
                        end   
                    end 
                elseif trialType == 2 || trialType == 4 
                    for z = 1:size(BdataToPlot{posCell},1)    
                         if isempty(BdataToPlot{posCell}{z,trialType}) == 0
                            for trial = 1:size(BdataToPlot{posCell}{z,trialType},2)      
                                if isempty(BdataToPlot{posCell}{z,trialType}{trial}) == 0 
                                    resamped_BdataToPlot{posCell}{z,trialType}{trial} = resample(BdataToPlot{posCell}{z,trialType}{trial},d2_len{2},d1_len{trialType});                                                                    
                                end 
                            end 
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
    
    
    for posCell = 1:length(BdataToPlot2)
        if isempty(BdataToPlot2{posCell}) == 0             
            for trialType = 1:size(BdataToPlot2{posCell},2)               
                for z = 1:size(BdataToPlot2{posCell},1)     
                     if isempty(BdataToPlot{posCell}{z,trialType}) == 0
                        for trial = 1:size(BdataToPlot2{posCell}{z,trialType},2) 
                            if isempty(BdataToPlot2{posCell}{z,trialType}{trial}) == 0 
                            resamped_BdataToPlot2{posCell}{z,trialType}{trial} = resample(BdataToPlot2{posCell}{z,trialType}{trial},d1_len{trialType},d2_len{trialType});                                                  
                            end 
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