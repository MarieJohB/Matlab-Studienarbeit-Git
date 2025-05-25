function cleaned_str = removeBracketsFromInput(input_str)
    % REMOVEBRACKETSFROMINPUT - Removes square brackets from input string
    %
    % This function removes square brackets [ ] from matrix input strings
    % to make them compatible with the str2matrix function
    %
    % Input:
    %   input_str - String with potential square brackets
    %
    % Output:
    %   cleaned_str - String with square brackets removed
    
    if isempty(input_str)
        cleaned_str = input_str;
        return;
    end
    
    % Remove all square brackets
    cleaned_str = strrep(input_str, '[', '');
    cleaned_str = strrep(cleaned_str, ']', '');
    
    % Trim any extra whitespace
    cleaned_str = strtrim(cleaned_str);
end