function [stackOut,BG_ROIboundData] = backgroundSubtraction(reg__Stacks)

disp('Now select a background ROI for background subtraction')
BG_ROIboundData = cell(1,length(reg__Stacks));
BG_ROIstacks = cell(1,length(reg__Stacks));
for stack = 1:length(reg__Stacks)
    data = reg__Stacks{stack};
    if stack == 1
        [ROI_stacks,xmins,ymins,widths,heights] = firstTimeCreateROIs(1,data);
        BG_ROIboundData{stack}{1} = xmins;
        BG_ROIboundData{stack}{2} = ymins;
        BG_ROIboundData{stack}{3} = widths;
        BG_ROIboundData{stack}{4} = heights;
        BG_ROIstacks{stack} = ROI_stacks;

    elseif stack > 1 
        BG_ROIboundData{stack}{1} = BG_ROIboundData{1}{1};
        BG_ROIboundData{stack}{2} = BG_ROIboundData{1}{2};
        BG_ROIboundData{stack}{3} = BG_ROIboundData{1}{3};
        BG_ROIboundData{stack}{4} = BG_ROIboundData{1}{4};

        xmins = BG_ROIboundData{stack}{1};
        ymins = BG_ROIboundData{stack}{2};
        widths = BG_ROIboundData{stack}{3};
        heights = BG_ROIboundData{stack}{4};
        [ROI_stacks] = make_ROIs_notfirst_time(data,xmins,ymins,widths,heights);
        BG_ROIstacks{stack} = ROI_stacks;
    end 
end 

% determine average pixel intensity of each frame in the control ROI
BGpixInt = cell(1,length(reg__Stacks));
for stack = 1:length(reg__Stacks)
    BGpixInt{stack} = mean(mean(BG_ROIstacks{stack}{1}));
end 

% do background subtraction 
stackOut = cell(1,length(reg__Stacks));
for stack = 1:length(reg__Stacks) 
    for frame = 1:size(BGpixInt{1},3)
        stackOut{stack}(:,:,frame) = (reg__Stacks{stack}(:,:,frame)-BGpixInt{stack}(:,:,frame));
    end 
end 

end 