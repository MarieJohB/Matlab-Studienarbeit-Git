function f = convert_string_to_function_handle(funcStr)
% CONVERT_STRING_TO_FUNCTION_HANDLE - Convert a string to a function handle
% This function converts a string expression to a MATLAB function handle,
% handling both time-dependent functions and constants.
%
% Parameters:
%   funcStr - String representation of a function or constant
%
% Returns:
%   f - Function handle representing the input string

% Check if the string is empty
if isempty(strtrim(funcStr))
    error('Function expression is empty.');
end

try
    % Check if the string contains 't' (time-dependent)
    if contains(funcStr, 't')
        % Create and test function handle
        f = str2func(['@(t)', funcStr]);
        test_result = f(0); % Test with a dummy value
    else
        % Try to convert to a number or evaluate as a MATLAB expression
        try
            % First try as a number
            eval_value = str2double(funcStr);
            if isnan(eval_value)
                % Then try as a MATLAB expression
                eval_value = eval(funcStr);
            end
            % Create function handle for constant value (adding 0*t to make it time-dependent)
            f = str2func(['@(t)', num2str(eval_value), ' + 0*t']);
        catch
            % If all else fails, try directly as a constant function
            f = str2func(['@(t)', funcStr, ' + 0*t']);
            % Test the function
            f(0);
        end
    end
catch ME
    error('Invalid function expression: %s', ME.message);
end
end