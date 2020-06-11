
% parse hdf file      
[parsedStruct,~]=csParser2;
teensyState = parsedStruct.teensyStates; 
sessionTime = parsedStruct.sessionTime;
vel_wheel_data=parsedStruct.velocity;
LED1_amp=parsedStruct.c1_amp;
LED2_amp=parsedStruct.c2_amp;
LED3_amp=parsedStruct.c3_amp;

[TrialTypes] = determineTrialType(LED1_amp,LED2_amp,LED3_amp);

state = 8;
%Find state start and end times 
[state_start_t,state_end_t] = find_state_bounds(teensyState,sessionTime,state);
%HACKY BUT IT WORKS 
%state_start_t = state_start_t - 5;

if length(state_start_t) > length(state_end_t)
    state_start_t(:,length(state_end_t)+1) = []; 
end 

state_start_t = floor(state_start_t);
state_end_t = floor(state_end_t);
state_start_f = (state_start_t)*(40000/sessionTime(end));
state_end_f = (state_end_t)*(40000/sessionTime(end));
statef_length = state_end_f-state_start_f; 
state_length = state_end_t-state_start_t;

state_start_t = state_start_t';
state_end_t = state_end_t';
state_length = state_length'; 
state_start_f = state_start_f';
state_end_f = state_end_f';
statef_length = statef_length'; 

HDFchart(:,1:2) = TrialTypes; 
HDFchart(:,3) = LED1_amp; 
HDFchart(:,4) = LED2_amp;
HDFchart(length(LED2_amp)+1-(length(LED2_amp)-length(state_start_t)):end,:) = [];
HDFchart(:,5) = state_start_t; 
HDFchart(:,6) = state_end_t; 
HDFchart(:,7) = state_length;
HDFchart(:,8) = state_start_f; 
HDFchart(:,9) = state_end_f; 
HDFchart(:,10) = statef_length;
%col_titles = ["Trials","Trial Type","LED2 Amp","LED3 Amp","Stim Start Time","Stim End Time","Stim Time Length","Stim Start Frame","Stim End Frame", "Stim Num Frames"];  
col_titles = ["Trials","Trial Type","LED1 Amp","LED2 Amp","Stim Start Time","Stim End Time","Stim Time Length","Stim Start Frame","Stim End Frame", "Stim Num Frames"];  
HDFchart = vertcat(col_titles,HDFchart); 


ResampedVel_wheel_data = resample(vel_wheel_data,40000,length(vel_wheel_data)); 

%get median wheel value 
WdataMed = median(ResampedVel_wheel_data);     
%compute Dv/v using means  
DVOV = (ResampedVel_wheel_data-WdataMed)./WdataMed;      
%get sliding baseline 
[WdataSlidingBL]=slidingBaseline(DVOV,floor(((40000/sessionTime(end)))*10),0.5); %0.5 quantile thresh = the median value                                   
%subtract sliding baseline from Dv/v
DVOVsubSBLs = DVOV-WdataSlidingBL;
%z-score wheel data 
zWData = zscore(DVOVsubSBLs);

stimTimes = [2;20];
stimTypeNum = 2;
[uniqueTrialData,uniqueTrialDataOcurr,indices,state_start_f,uniqueTrialDataTemplate] = separateTrialTypes(TrialTypes,state_start_f,state_end_f,stimTimes,1,(40000/sessionTime(end)),stimTypeNum); 

%sort wheel data                    
dataParseType = 0;
sec_before_stim_start = 20;
sec_after_stim_end = 20;
FPstack = (40000/sessionTime(end));
Tdata = zWData;

sortedTrials = cell(1,length(uniqueTrialDataOcurr));
for indGroup = 1:length(uniqueTrialDataOcurr)     
    stimOnFrames = uniqueTrialData(indGroup,2);

    if dataParseType == 0    
        counter = 1; 
        trial_start = zeros(1,counter);
        trial_end = zeros(1,counter);
        for ind = 1:length(indices{indGroup})
            %determine trial start and end frames
            trial_start(counter) = floor(state_start_f(indices{indGroup}(ind)) - (sec_before_stim_start*FPstack));
            trial_end(counter) = trial_start(counter) + ((stimOnFrames+((sec_after_stim_end+sec_before_stim_start)*FPstack)));            
            trial_end = trial_end';             
            counter = counter + 1;
        end     

        %sort data 
        for ind = 1:length(indices{indGroup})       
            if trial_start(ind) > 0 && trial_end(ind) < length(Tdata) 
                sortedTrials{indGroup}{ind} = Tdata(trial_start(ind):trial_end(ind));
            elseif trial_start(ind) < 0 
                indices{indGroup}(ind) = [];
            end 
        end   

    elseif dataParseType == 1           
        counter = 1;            
            for ind = 1:length(indices{indGroup})
                %determine trial start and end frames
                trial_start(counter) = floor(state_start_f(indices{indGroup}(ind)) - (sec_before_stim_start*FPstack));
                trial_end(counter) = trial_start(counter) + ((stimOnFrames+((sec_after_stim_end+sec_before_stim_start)*FPstack)));            
                trial_end = trial_end';               
                counter = counter + 1;
            end
    end

    %sort data 
    for ind = 1:length(indices{indGroup})                   
        if trial_start(ind) > 0 && trial_end(ind) < length(Tdata) 
            sortedWheelData{indGroup}{ind} = Tdata(trial_start(ind):trial_end(ind));
        end  
    end 
end 


%sort wheel data into correct spot based on trial type 
for trialType = 1:4
    if ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows') == 1
        [~, idxStart] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialData,'rows');  
        [~, idxFin] = ismember(uniqueTrialDataTemplate(trialType,:),uniqueTrialDataTemplate,'rows');
                  
        wheelDataToPlot{idxFin} = sortedWheelData{idxStart};                    
    end     
end

FPS = (40000/sessionTime(end));
numZplanes = 1;