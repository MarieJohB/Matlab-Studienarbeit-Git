function user_function = get_function_input()

target_value = [];

while isempty(target_value)
        % request input
        prompt = {'Enter the target value (use "t" for time-dependent values, e.g., "5" or "sin(t)"): '};
        dlgtitle = 'Target Value Input';
        dims = [1 50];
        answer = inputdlg(prompt, dlgtitle, dims);
        
        % check if user cancelled 
        if isempty(answer)
            disp('Operation cancelled by user.');
            return;
        end
        
        % Validate the entered input
        try
            % replacing commas with periods
            target_value_str = strrep(answer{1}, ',', '.');
            
            % check whether string is empty
            if isempty(strtrim(target_value_str))
                error('Invalid input. Please enter a valid expression.');
            end
            
            % does input contain 't' --> time dependent function
            if contains(target_value_str, 't')
                % create and test function
                test_func = str2func(['@(t)', target_value_str]);
                test_result = test_func(0); % dummy value for test
                % if test does not fail: save function in array
                target_value = test_func;
            else
                % no 't' contained --> constant value / not time dependent 
                eval_value = str2double(target_value_str);
                if isnan(eval_value)
                    error('Invalid input. Please enter a valid expression.');
                end
                % If it's a constant value, create a function that returns the constant + 0*t
                target_value = str2func(['@(t)', num2str(eval_value), ' + 0*t']);
            end
        catch
            uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
            target_value = [];
        end
end

% show function
disp('Stored functions:');
disp(target_value);

% return the created function
user_function = target_value;

end