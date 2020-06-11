function [DiffIms] = makeDiffIm(inputArray)

DiffIms = cell(1,size(inputArray,1));
for rows = 1:size(inputArray,1)
    if size(inputArray,1) == 1
        for cols = 1:size(inputArray{rows},2)
            DiffIms{rows}{cols} = inputArray{rows}{cols}(:,:,end) - inputArray{rows}{cols}(:,:,1); 
        end 
    elseif size(inputArray,1) > 1
        DiffIms{rows} = inputArray{rows}(:,:,end) - inputArray{rows}(:,:,1); 

    end 

end 



end 