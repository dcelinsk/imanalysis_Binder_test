function [subStacks] = splitStacks(imageStack,splitType,numSplits)
subStackLength = size(imageStack,3)/numSplits; 

if splitType == 0 
    subStacks = cell(1,numSplits);
    for subStack = 1:numSplits
        for subFrame = 1:subStackLength
            if subStack == 1 
                subStacks{subStack}(:,:,subFrame) = imageStack(:,:,subFrame);
            elseif subStack > 1 
                frame = subFrame + subStackLength*(subStack-1);
                subStacks{subStack}(:,:,subFrame) = imageStack(:,:,frame);
            end 
        end 
    end 
               
elseif splitType == 1   
    subStacks = cell(1,numSplits);
    for subStack = 1:numSplits
        for subFrame = 1:subStackLength
            if subFrame == 1 
                frame = subStack; 
            end 
            subStacks{subStack}(:,:,subFrame) = imageStack(:,:,frame);
            frame = frame + numSplits;
        end 
    end 
end 

end 