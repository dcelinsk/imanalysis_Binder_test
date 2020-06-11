%% set paramaters 
clearvars 

animal = 'SF56_20190718';
ROI = 'ROI2';
videoNum = 1;

G1orR2 = 1;
if G1orR2 == 1
    col = 3;
    color = 'green';
elseif G1orR2 == 2
    col = 4; 
    color = 'red';
end 

%% get the data
temp = matfile(sprintf('%s_%s_%d_regIms_%s.mat',animal,ROI,videoNum,color));
regStacks = temp.regStacks; 



%% average registered imaging data across planes in Z 
inputStackArray = zeros(size(regStacks{2,col}{1},1),size(regStacks{2,col}{1},2),size(regStacks{2,col}{1},3),size(regStacks{2,col},2));
for Z = 1:size(regStacks{2,col},2)
    inputStackArray(:,:,:,Z) = regStacks{2,col}{Z};
end 
stack1 = mean(inputStackArray,4);

%% remove background noise 
stack = zeros(size(regStacks{2,col}{1},1),size(regStacks{2,col}{1},2),size(regStacks{2,col}{1},3));
for frame = 1:size(stack1,3)
    stack(:,:,frame) = stack1(:,:,frame) - mean(mean(stack1(:,:,frame)));
end 


%% save the stack 
% dir = sprintf('D:/70kD_RhoB/DAT-Chrimson-GCaMP/%s/%s_%s_%d_%sIms_avAcrossZ_BGnoiseNotRemoved',animal,animal,ROI,videoNum,color);
% mkdir(dir);
% 
% for frame = 1:size(stack,3)
%     folder = sprintf('%s/%sIm_%d.tif',dir,color,frame');
%     imwrite(mat2gray(stack(:,:,frame)),folder);
% end 

%% create cumulative pixel intensity stack
% 
% for frame = 1:size(stack2,3)
%     if frame == 1 
%         stack3(:,:,frame) = stack2(:,:,frame);
%     elseif frame > 1 
%         stack3(:,:,frame) = stack2(:,:,frame)+stack3(:,:,frame-1);
%     end 
% end 









