function [dffStacks,CumDffStacks,CumStacks] = makeCumPixStacksPerTrial(sortedStacks,FPS,numZplanes,sec_before_stim_start)

%make baseline image 
TtypeBaselineIms = cell(1,length(sortedStacks));
for zStack = 1:length(sortedStacks)
    for trialType = 1:length(sortedStacks{zStack})
        for stack = 1:length(sortedStacks{zStack}{trialType})
            TtypeBaselineIms{zStack}{trialType}{stack} = mean(sortedStacks{zStack}{trialType}{stack}(:,:,1:ceil((FPS/numZplanes)*(sec_before_stim_start/2))),3);
        end   
    end 
end 

%create df/f stack relative to the first 20 sec of the exp (before first
%stim) 
dffStacks = cell(1,length(sortedStacks));
for zStack = 1:length(sortedStacks)
    for trialType = 1:length(sortedStacks{zStack})
        for stack = 1:length(sortedStacks{zStack}{trialType})
            dffStacks{zStack}{trialType}{stack} = (double(sortedStacks{zStack}{trialType}{stack})-repmat(TtypeBaselineIms{zStack}{trialType}{stack},[1 1 size(double(sortedStacks{zStack}{trialType}{stack}),3)]))./repmat(TtypeBaselineIms{zStack}{trialType}{stack},[1 1 size(double(sortedStacks{zStack}{trialType}{stack}),3)]);
            dffStacks{zStack}{trialType}{stack} = int16(dffStacks{zStack}{trialType}{stack});
        end   
    end 
end 

%create cumulative df/f stacks normalized to baseline
CumDffStacks = cell(1,length(sortedStacks));
for zStack = 1:length(sortedStacks)
    for trialType = 1:length(sortedStacks{zStack})
        for stack = 1:length(sortedStacks{zStack}{trialType})
            for frame = 1:size(sortedStacks{zStack}{trialType}{stack},3)    
                if frame == 1
                    CumDffStacks{zStack}{trialType}{stack}(:,:,frame) = dffStacks{zStack}{trialType}{stack}(:,:,frame);
                elseif frame > 1
                    CumDffStacks{zStack}{trialType}{stack}(:,:,frame) = CumDffStacks{zStack}{trialType}{stack}(:,:,frame-1) + dffStacks{zStack}{trialType}{stack}(:,:,frame); %sum(dffStacksFirst20s{stack}(:,:,1:frame),3);
                end   
            end 
        end 
    end 
end 

%create cumulative pixel intensity stacks
CumStacks = cell(1,length(sortedStacks));
for zStack = 1:length(sortedStacks)
    for trialType = 1:length(sortedStacks{zStack})
        for stack = 1:length(sortedStacks{zStack}{trialType})
            for frame = 1:size(sortedStacks{zStack}{trialType}{stack},3) 

                if frame == 1
                    CumStacks{zStack}{trialType}{stack}(:,:,frame) = sortedStacks{zStack}{trialType}{stack}(:,:,frame);
                elseif frame > 1
                    CumStacks{zStack}{trialType}{stack}(:,:,frame) = CumStacks{zStack}{trialType}{stack}(:,:,frame-1) + sortedStacks{zStack}{trialType}{stack}(:,:,frame); %sum(dffStacksFirst20s{stack}(:,:,1:frame),3);
                end 
            end 
        end 
    end                 
end 