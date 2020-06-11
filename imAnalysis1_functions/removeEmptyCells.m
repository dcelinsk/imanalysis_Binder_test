function [sortedStacks,indices,emptyTrialTypes] = removeEmptyCells(sortedStacks,indices)
emptyTrialTypes = cell(1,length(sortedStacks));
emptyTrials = cell(1,length(sortedStacks));
for zStack = 1:length(sortedStacks)
    %find empty trialTypes 
    emptyTrialTypes{zStack} = cellfun(@isempty, sortedStacks{zStack});   
    for trialType = 1:length(sortedStacks{zStack}) 
        if emptyTrialTypes{zStack}(trialType) == 0           
            %find empty trials 
            emptyTrials{zStack}{trialType} = cellfun(@isempty, sortedStacks{zStack}{trialType});
            %remove empty trials
            sortedStacks{zStack}{trialType} = sortedStacks{zStack}{trialType}(~cellfun('isempty',sortedStacks{zStack}{trialType}));    
            for trial = 1:length(sortedStacks{zStack}{trialType})
                if zStack == 1 && emptyTrials{zStack}{trialType}(trial) == 1 
                    indices{trialType}(trial) = []; 
                end 
            end        
        end 
    end 
end
end 