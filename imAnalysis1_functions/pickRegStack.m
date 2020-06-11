function [reg_Stacks] = pickRegStack(regStacks,regTypeDim,regTypeTemp)
if regTypeDim == 0 
    if regTypeTemp == 0 
        reg_Stacks = regStacks{2,2};
    elseif regTypeTemp == 1 
        reg_Stacks = regStacks{2,1};
    end 
elseif regTypeDim == 1 
    if regTypeTemp == 0 
        reg_Stacks = regStacks{2,4};
    elseif regTypeTemp == 1 
        reg_Stacks = regStacks{2,3};
    end 
end 