function [subStacks] = reorgVolStack(imageStack,splitType,numSplits)
subStackLength = floor(size(imageStack,3)/numSplits); 

if splitType == 0 
    subStacks = zeros(size(imageStack,1),size(imageStack,2),numSplits,subStackLength);
    for subStack = 1:numSplits
        for subFrame = 1:subStackLength
            if subStack == 1 
                subStacks(:,:,subStack,subFrame) = imageStack(:,:,subFrame);
            elseif subStack > 1 
                frame = subFrame + subStackLength*(subStack-1);
                subStacks(:,:,subStack,subFrame) = imageStack(:,:,frame);
            end 
        end 
    end 
               
elseif splitType == 1   
    subStacks = zeros(size(imageStack,1),size(imageStack,2),numSplits,subStackLength);
    for subStack = 1:numSplits
        for subFrame = 1:subStackLength
            if subFrame == 1 
                frame = subStack; 
            end 
            subStacks(:,:,subStack,subFrame) = imageStack(:,:,frame);
            frame = frame + numSplits;
        end 
    end 
end 

end 