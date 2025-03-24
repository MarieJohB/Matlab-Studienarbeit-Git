function str = conditionalString(condition, trueStr, falseStr)
    % Helper function to return a string based on a condition
    if condition
        str = trueStr;
    else
        str = falseStr;
    end
end