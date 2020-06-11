function [dataToPlot] = catData3(dataToPlot,dataToPlot2)


for z = 1:length(dataToPlot)
    for trialType = 1:size(dataToPlot{z},2)  
        startSize = size(dataToPlot{z}{trialType},2);
        if trialType <= size(dataToPlot2{z},2) 
            if isempty(dataToPlot2{z}{trialType}) == 0 
                for trial = 1:size(dataToPlot2{z}{trialType},2)   
                    %dataToPlot{z}{trialType}{size(dataToPlot{z}{trialType},2)+1:size(dataToPlot2{z}{trialType},2)+size(dataToPlot{z}{trialType},2),:} = dataToPlot2{z}{trialType};  
                    %dataToPlot3{z}{trialType} = dataToPlot2{z}{trialType};

                    dataToPlot{z}{trialType}{startSize+trial} = dataToPlot2{z}{trialType}{trial}; 


                end 

            end
        elseif trialType > size(dataToPlot2{z},2) 
            
            dataToPlot{z}{trialType} = dataToPlot{z}{trialType}; 
            
        end 
    end 
  
end 



end 