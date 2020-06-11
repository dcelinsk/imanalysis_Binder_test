
function [ROI_stacks,xmins,ymins,widths,heights] = firstTimeCreateROIs(num_ROIs,data)
    
    ROI_stacks = cell(1,num_ROIs);
    xmins = zeros(1,num_ROIs);
    ymins = zeros(1,num_ROIs);
    widths = zeros(1,num_ROIs);
    heights = zeros(1,num_ROIs);
    for x = 1:num_ROIs
        [ROI_stack,xmin,ymin,width,height] = make_ROIs_first_time(data);
        ROI_stacks{x} = ROI_stack; 
        xmins(x) = xmin;
        ymins(x) = ymin;
        widths(x) = width;
        heights(x) = height; 
    end 
end 