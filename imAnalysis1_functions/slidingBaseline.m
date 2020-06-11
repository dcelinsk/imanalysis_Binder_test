function [out]=slidingBaseline(data,windowSize,quantileThresh)


%pre-allocate the vector and then map the data into the pieces between the
%window fragments, then map the first window/2 values into the first pad
%and the last window/2 values into the last pad.
out=zeros(size(data));


if mod(windowSize,2)
    windowSize=windowSize+1;
else
end

halfInd=fix(windowSize/2);



padD=zeros(size(data,1),size(data,2)+windowSize);
%disp(size(padD))

padD(:,halfInd+1:end-halfInd)=data;
padD(:,1:halfInd)=data(:,halfInd+1:(halfInd+1)+(halfInd-1));
padD(:,end-(halfInd-1):end)=data(:,end-(halfInd-1):end);
%disp(size(padD))

startQuant=halfInd;
if size(data,1)>1
    parfor n=1:size(data,2)
        td=padD(:,(n+startQuant)-startQuant:(n+startQuant)+(startQuant-1));
        aa=quantile(td,quantileThresh,2);
        out(:,n)=diag(aa);
    end
else
    parfor n=1:size(data,2)
        td=padD(:,(n+startQuant)-startQuant:(n+startQuant)+(startQuant-1));
        aa=quantile(td,quantileThresh);
        out(:,n)=aa;
    end
end


