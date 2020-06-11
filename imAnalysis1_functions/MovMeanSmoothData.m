function [filtData] = MovMeanSmoothData(data,filtTime,FPS)

filter_rate = FPS*filtTime; 

filtData = smoothdata(data,2,'movmean',filter_rate);

end 

