function [Qanswer] = getUserInput(userInput,Qstring)
    %find the question
    Q = strfind(userInput,Qstring);
    isOne = cellfun(@(x)isequal(x,1),Q);
    [row,col] = find(isOne); 
    
    if Qstring == "Where Photos Are Found" || Qstring == "imAnalysis1_functions Directory" || Qstring == "imAnalysis2_functions Directory"
        Qanswer = userInput(row,col+1);  
    elseif Qstring ~= "Where Photos Are Found"    
        %check if there are spaces in answer (as w/ the angles)
        if sum(isspace(userInput(row,col+1))) > 0 
            %get angle answer
            Qanswer = split(userInput(row,col+1));
            Qanswer = str2double(Qanswer);

        elseif sum(isspace(userInput(row,col+1))) == 0
            %get the answer 
            Qanswer = str2double(userInput(row,col+1));
        end 
    end 

end 