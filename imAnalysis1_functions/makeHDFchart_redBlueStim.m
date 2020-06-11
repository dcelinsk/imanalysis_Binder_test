function [HDFchart,state_start_f,state_end_f,FPS,vel_wheel_data,TrialTypes] = makeHDFchart_redBlueStim(state,framePeriod)

% parse hdf file      
[parsedStruct,~]=csParser2;
teensyState = parsedStruct.teensyStates; 
sessionTime = parsedStruct.sessionTime;
vel_wheel_data=parsedStruct.velocity;
LED1_amp=parsedStruct.c1_amp;
LED2_amp=parsedStruct.c2_amp;
LED3_amp=parsedStruct.c3_amp;

[TrialTypes] = determineTrialType(LED1_amp,LED2_amp,LED3_amp);

%Find state start and end times 
[state_start_t,state_end_t] = find_state_bounds(teensyState,sessionTime,state);
%HACKY BUT IT WORKS 
%state_start_t = state_start_t - 5;

if length(state_start_t) > length(state_end_t)
    state_start_t(:,length(state_end_t)+1) = []; 
end 

FPS = 1/framePeriod; 

state_start_t = floor(state_start_t);
state_end_t = floor(state_end_t);
state_start_f = (state_start_t)*FPS;
state_end_f = (state_end_t)*FPS;
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
end 