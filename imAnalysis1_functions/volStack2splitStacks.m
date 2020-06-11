function [VZstacks] = volStack2splitStacks(regVolStack)

VZstacks = cell(1,size(regVolStack,3));
for subStack = 1:size(regVolStack,3)
    for subFrame = 1:size(regVolStack,4)
        VZstacks{subStack}(:,:,subFrame) = regVolStack(:,:,subStack,subFrame);
    end 
end 
       

end 