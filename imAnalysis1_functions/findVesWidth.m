function [vessel_diam,boundaries,area] = findVesWidth(BWstacks)

vessel_diam = zeros(1,size(BWstacks,3));
area = zeros(1,size(BWstacks,3));
boundaries = cell(1,size(BWstacks,3));
for y = 1:size(BWstacks,3)
        %determine number of white pixels per row of mask (BW)
        white_pix = flipud(sum(BWstacks(:,:,y) == 1, 2));
        % determine the area by summing the counts of all white pixels
        % across each row of an image
        area(y) = sum(white_pix); % Dmitrijs
        %find mean number of white pixels per row
        av_white_pix = mean(white_pix);
        %find diameter of largest region
        vessel_diam(y) = av_white_pix;
        %make outline
        [B,~] = bwboundaries(BWstacks(:,:,y),'noholes');
        boundaries{y} = B;
end 


end 