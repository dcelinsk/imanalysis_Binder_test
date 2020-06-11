function [uniqueTrialData,uniqueTrialDataOcurr,indices,state_start_f,uniqueTrialDataTemplate] = separateTrialTypes(TrialTypes,state_start_f,state_end_f,stimTimes,numZplanes,FPS,stimTypeNum)

Ttypes = TrialTypes(:,2);
%determine theoretical trial lengths (in frames) 
theoTLengths = zeros(1,length(stimTimes));
for L = 1:length(stimTimes)
    theoTLengths(L) = (FPS/numZplanes)*stimTimes(L);
end 
theoTLengths = floor(theoTLengths); 

%determine the actual trial lengths (in frames) = col 1. col 2 = num of
%occurances 
trialLengths = state_end_f - state_start_f; 
lengths = unique(trialLengths); 
lengths(:,2) = histc(trialLengths(:),lengths);    
%identify trial lengths that should be the same (kmeans clustering)
if size(lengths,1) > 1 
    [lengthGroups,centroids] = kmeans(lengths(:,1),2);
elseif size(lengths,1) == 1
    lengthGroups = 1;
    centroids = trialLengths;
end 

UniqueLengthGroups = unique(lengthGroups);

%% remove trials that are too long (more than 10% off) = where there was probably a mechanical/triggering error 
groupRows = cell(1,length(UniqueLengthGroups));
diffArray = zeros(length(UniqueLengthGroups),length(theoTLengths));
 for uniqueGroup = 1:length(UniqueLengthGroups)
    [groupRow, ~] = find(lengthGroups == UniqueLengthGroups(uniqueGroup));
    groupRows{uniqueGroup} = groupRow;
    for numLengths = 1:length(theoTLengths)
        diffArray(uniqueGroup,numLengths) = abs(centroids(uniqueGroup) - theoTLengths(numLengths));
    end 

    [~, minCol] = min(diffArray(uniqueGroup,:));
    centroids(uniqueGroup,2) = theoTLengths(minCol);

     for ind = 1:length(groupRows{uniqueGroup})
        if lengths(groupRows{uniqueGroup}(ind)) > (0.1*centroids(uniqueGroup,1))+centroids(uniqueGroup,1)

            state_start_f(find(trialLengths == lengths(groupRows{uniqueGroup}(ind)))) = [];
            state_end_f(find(trialLengths == lengths(groupRows{uniqueGroup}(ind)))) = [];          
            Ttypes(find(trialLengths == lengths(groupRows{uniqueGroup}(ind))),:) = [];
        end 
     end 
 end 
%%
%determine the actual trial lengths (in frames) again 
state_start_f = floor(state_start_f); 
state_end_f = ceil(state_end_f);
trialLengths = state_end_f - state_start_f; 
lengths = unique(trialLengths); 
lengths(:,2) = histc(trialLengths(:),lengths);
%identify trial lengths that should be the same (kmeans clustering)
if size(lengths,1) > 1 
    [lengthGroups,centroids] = kmeans(lengths(:,1),2);
elseif size(lengths,1) == 1
    lengthGroups = 1;
end 

%standardize length of trials 
for uniqueGroup = 1:length(UniqueLengthGroups)
    [groupRow, ~] = find(lengthGroups == UniqueLengthGroups(uniqueGroup));
    groupRows{uniqueGroup} = groupRow;
    for numLengths = 1:length(theoTLengths)
        diffArray(uniqueGroup,numLengths) = abs(centroids(uniqueGroup) - theoTLengths(numLengths));
    end 

    [~, minCol] = min(diffArray(uniqueGroup,:));
    centroids(uniqueGroup,2) = theoTLengths(minCol);

    for ind = 1:length(groupRows{uniqueGroup})
        if lengths(groupRows{uniqueGroup}(ind)) ~= centroids(uniqueGroup,2)
            state_end_f(find(trialLengths == lengths(groupRows{uniqueGroup}(ind)))) = state_start_f(find(trialLengths == lengths(groupRows{uniqueGroup}(ind)))) + centroids(uniqueGroup,2);               
        end 
    end 
end 

%determine the actual trial lengths (in frames) again 
state_start_f = floor(state_start_f); 
state_end_f = ceil(state_end_f);
trialLengths = state_end_f - state_start_f; 
lengths = unique(trialLengths); 
lengths(:,2) = histc(trialLengths(:),lengths);

if length(Ttypes) ~= length(trialLengths)
    if length(Ttypes) < length(trialLengths)
        trialLengths = trialLengths(1:length(Ttypes));
    elseif length(Ttypes) > length(trialLengths)
        Ttypes = Ttypes(1:length(trialLengths));
    end 
end 
trialData = horzcat(Ttypes,trialLengths);    
[uniqueTrialData,~,ib] = unique(trialData,'rows');
uniqueTrialDataOcurr = accumarray(ib, 1);
indices = accumarray(ib, find(ib), [], @(rows){rows});  %the find(ib) simply generates (1:size(a,1))'

%make a template for sorting trial types later on 
stimTypes = (1:stimTypeNum);
stimTypesU = repmat(stimTypes,length(theoTLengths),1);
stimTypesU = stimTypesU(:);
theoTLengths = theoTLengths';
theoTLengthsU = repmat(theoTLengths,stimTypeNum,1);

uniqueTrialDataTemplate(:,1) = stimTypesU;
uniqueTrialDataTemplate(:,2) = theoTLengthsU;

end 