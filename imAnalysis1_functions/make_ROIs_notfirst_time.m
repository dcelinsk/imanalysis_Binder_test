function [ROI_stacks] = make_ROIs_notfirst_time(segmented_data,xmins,ymins,widths,heights)

ROI_stacks = cell(1,length(xmins));
for x = 1: length(xmins)
    for i=1: size(segmented_data,3) 
        %You can now go ahead and crop the image:
        ROI_stacks{x}(:,:,i) = imcrop(segmented_data(:,:,i), [xmins(x) ymins(x) widths(x) heights(x)]);
    end    
end 
end 