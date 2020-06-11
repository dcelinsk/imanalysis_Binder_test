%% get just the data you need 
temp = matfile('SF58_20190627_70FITC_ROI2_V3-5_BBB.mat');
SF58_ROI2_Bdata = temp.BdataToPlot; %EDIT 
FPS = temp.FPS;
if FPS > 20 
    FPS = FPS / 3 ; 
end 
SF58_ROI2_FPS = FPS; %EDIT
clear FPS temp

%% organize data into giant cell array - iterate through mice, ROI, Z, trialType, trial, FOV 
mouse = 4;
ROI = 2; 
data{mouse}{ROI} = SF58_ROI2_Bdata;
FPSdata(mouse,ROI) = SF58_ROI2_FPS;

%% upsample data 
maxFPS = max(max(FPSdata));
Rdata = cell(1,length(data));
for mouse = 1:length(data)
    for ROI = 1:length(data{mouse})
        for Z = 1:length(data{mouse}{ROI})
            for trialType = 1:length(data{mouse}{ROI}{Z})
                for trial = 1:length(data{mouse}{ROI}{Z}{trialType})
                    for FOV = 1:length(data{mouse}{ROI}{Z}{trialType}{trial})
                        if trialType == 1 || trialType == 3
                            goalLen = round(maxFPS*42);
                            Rdata{mouse}{ROI}{Z}{trialType}{trial}{FOV} = resample(data{mouse}{ROI}{Z}{trialType}{trial}{FOV},goalLen,length(data{mouse}{ROI}{Z}{trialType}{trial}{FOV}));
                        elseif trialType == 2 || trialType == 4
                            goalLen = round(maxFPS*60);
                            Rdata{mouse}{ROI}{Z}{trialType}{trial}{FOV} = resample(data{mouse}{ROI}{Z}{trialType}{trial}{FOV},goalLen,length(data{mouse}{ROI}{Z}{trialType}{trial}{FOV}));
                        end 
                    end 
                end 
            end 
        end 
    end 
end 

%% create cell array of timeseries objects 

% TSdata = cell(1,length(data));
% for mouse = 1:length(data)
%     for ROI = 1:length(data{mouse})
%         for Z = 1:length(data{mouse}{ROI})
%             for trialType = 1:length(data{mouse}{ROI}{Z})
%                 for trial = 1:length(data{mouse}{ROI}{Z}{trialType})
%                     for FOV = 1:length(data{mouse}{ROI}{Z}{trialType}{trial})
%                         Tval = ((1/FPSdata{mouse}{ROI}):(1/FPSdata{mouse}{ROI}):(length(data{mouse}{ROI}{Z}{trialType}{trial}{FOV}))/FPSdata{mouse}{ROI});
%                         TSdata{mouse}{ROI}{Z}{trialType}{trial}{FOV} = timeseries(data{mouse}{ROI}{Z}{trialType}{trial}{FOV},Tval);
%                     end 
%                 end 
%             end 
%         end 
%     end 
% end 

%% clear unescessary values 
clearvars -except data FPSdata TSdata Rdata
