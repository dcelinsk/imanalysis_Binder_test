function [sortedTrials,indices] = eventTriggeredAverages2(Tdata,state_start_f,indices,uniqueTrialData,uniqueTrialDataOcurr,numZplanes)

    trialStart = floor(state_start_f/numZplanes);
  
    sortedTrials = cell(1,length(uniqueTrialDataOcurr));
    for indGroup = 1:length(uniqueTrialDataOcurr)    
        stimOnFrames = uniqueTrialData(indGroup,2);
            
        counter = 1;            
        for ind = 1:length(indices{indGroup})
            %determine trial start and end frames for stims 
            trial_start(counter) = floor(trialStart(indices{indGroup}(ind)));
            trial_end(counter) = floor(trial_start(counter) + stimOnFrames);            
             
            
            
            trialEnd(indices{indGroup}(ind)) = trial_end(counter);
            counter = counter + 1;
        end
        
%         trialEnd = trialEnd';
        baselineStart = 1;
        baselineStart = vertcat(baselineStart,(trialEnd+1)');
        baselineEnd = trialStart-1;
        
        
        %sort data 
        for ind = 1:length(indices{indGroup})                   
            if trial_start(ind) > 0 && trial_end(ind) < length(Tdata) 
                sortedTrials{indGroup}{ind} = Tdata(trial_start(ind):trial_end(ind));
            end  
        end 
    end 
    
    putHere = length(sortedTrials)+1;
    counter2 = 1;
    for basePer = 1:length(baselineEnd)
        if baselineStart(basePer) < length(Tdata) && baselineEnd(basePer) < length(Tdata)
            sortedTrials{putHere}{counter2} = Tdata(baselineStart(basePer):baselineEnd(basePer)); 
            counter2 = counter2 + 1;
        end 
    end 
end 
    