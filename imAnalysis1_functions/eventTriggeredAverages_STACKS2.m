function [sortedTrials,indices] = eventTriggeredAverages_STACKS2(Tdata,state_start_f,FPS,indices,uniqueTrialData,uniqueTrialDataOcurr,userInput,numZplanes)
    

    FPstack = FPS/numZplanes;
    [dataParseType] = getUserInput(userInput,"What data do you need? Peristimulus epoch = 0. Stimulus epoch = 1.");
    
   
    if dataParseType == 0 
        sec_before_stim_start = getUserInput(userInput,"How many seconds before the stimulus starts do you want to plot?");
        sec_after_stim_end = getUserInput(userInput,"How many seconds after stimulus end do you want to plot?");
    elseif dataParseType == 1 
        sec_before_stim_start = 0;
        sec_after_stim_end = 0; 
    end 
  
    sortedTrials = cell(1,length(uniqueTrialDataOcurr));
    for indGroup = 1:length(uniqueTrialDataOcurr)     
        stimOnFrames = uniqueTrialData(indGroup,2);
        
        if dataParseType == 0    
            counter = 1; 
            trial_start = zeros(1,1);
            trial_end = zeros(1,1);
            for ind = 1:length(indices{indGroup})
                %determine trial start and end frames
                trial_start(counter) = floor(state_start_f(indices{indGroup}(ind)) - (sec_before_stim_start*FPstack));
                trial_end(counter) = trial_start(counter) + ((stimOnFrames+((sec_after_stim_end+sec_before_stim_start)*FPstack)));            
                trial_end = trial_end';            
                counter = counter + 1;
            end     
            
            %sort data and remove irrelevant indices 
            for ind = 1:length(indices{indGroup})       
                if trial_start(ind) > 0 && trial_end(ind) < length(Tdata) 
                    sortedTrials{indGroup}{ind} = Tdata(:,:,trial_start(ind):trial_end(ind));
                end 
            end   
  
        elseif dataParseType == 1           
            counter = 1;            
                for ind = 1:length(indices{indGroup})
                    %determine trial start and end frames
                    trial_start(counter) = floor(state_start_f(indices{indGroup}(ind)));
                    trial_end(counter) = trial_start(counter) + stimOnFrames;            
                    trial_end = trial_end';             
                    counter = counter + 1;
                end
        end
        
        %sort data 
        for ind = 1:length(indices{indGroup})                   
            if trial_start(ind) > 0 && trial_end(ind) < length(Tdata) 
                sortedTrials{indGroup}{ind} = Tdata(:,:,trial_start(ind):trial_end(ind));
            end  
        end 
    end 
end 
    