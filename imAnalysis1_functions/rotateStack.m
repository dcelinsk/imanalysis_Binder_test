function [rotStacks,rotateImAngle] = rotateStack(reg_Stacks)

%% rotate image? 
alreadyRot = 0;
imshow(mean(reg_Stacks{1},3),[0,1800]);
rotateIm = input("Do you need to rotate the stack? Yes = 1. No = 0. ");    

while rotateIm == 1
    rotateImAngle = input("How much do you need to rotate the image? ");
    %rotate display image
    rotImage = imrotate(reg_Stacks{1}(:,:,1),rotateImAngle);
    %show the rotated display image
    imshow(rotImage,[0,1000]);    
    rotateImAsk = input("Does the image need to be rotated again? Yes = 1. No = 0. ");
    
    %if rotation is good, rotate all z stacks  
    if rotateImAsk == 0
        rotStacks = cell(1,length(reg_Stacks));
        for stack = 1:length(reg_Stacks)  
           %if the rotation is finally good, rotate the stack 
            for frame = 1:size(reg_Stacks{stack},3)
                rotStacks{stack}(:,:,frame) = imrotate(reg_Stacks{stack}(:,:,frame),rotateImAngle);
            end 
        end 
        alreadyRot = 1;
        rotateIm = 0;
        
    %if the rotation isn't good, start over 
    elseif rotateImAsk == 1
        rotateIm = 1; 
    end 
end 

if rotateIm == 0 && alreadyRot == 0
    rotStacks = reg_Stacks;
    rotateImAngle = 0; 
end 


end 