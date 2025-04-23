function formatted_str = formatNumericValues(str)
    % Format numeric values in a string to have at most 5 decimal places
    % This is a simplified approach that works with basic mathematical expressions
    
    % Find all numbers in the string (including decimal and scientific notation)
    [numbers, positions] = regexp(str, '-?\d+\.\d+|-?\d+\.?\d*e[-+]?\d+', 'match', 'start');
    
    % Process the string from end to beginning to avoid position shifts
    for i = length(numbers):-1:1
        num_str = numbers{i};
        num_val = str2double(num_str);
        
        % Round to 5 decimal places
        rounded_val = round(num_val * 10^5) / 10^5;
        new_num_str = num2str(rounded_val);
        
        % Replace in the original string
        str = [str(1:positions(i)-1), new_num_str, str(positions(i)+length(num_str):end)];
    end
    
    formatted_str = str;
end