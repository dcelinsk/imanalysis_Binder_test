function play_mask_over_roi_stack(stack,boundaries,framesToShow,maxPoint)


for frame = 1:framesToShow
    imshow(stack(:,:,frame),[0,maxPoint])
    hold on; 
    for bounds = 1:size(boundaries{frame},1)
        plot(boundaries{frame}{bounds,1}(:,2), boundaries{frame}{bounds,1}(:,1), 'r', 'LineWidth', 2)
    end 
end 
end 

