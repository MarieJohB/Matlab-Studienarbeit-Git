function target_value = get_target_value()
    % Initialize target_value to be empty
    target_value = [];
    
    % Loop until valid input is provided or user cancels the dialog
    while isempty(target_value)
        % Prompt for target value (constant or time-dependent)
        prompt = {'Enter the target value (use "t" for time-dependent values, e.g., "5" or "sin(t)"): '};
        dlgtitle = 'Target Value Input';
        dims = [1 50];
        answer = inputdlg(prompt, dlgtitle, dims);

        % Check if user cancelled the dialog
        if isempty(answer)
            disp('Operation cancelled by user.');
            return;
        end

        % Validate input
        try
            % Replace commas with periods
            target_value_str = strrep(answer{1}, ',', '.');

            % Test the input
            t = 0; % Dummy value for validation
            target_value = eval(target_value_str);
            if isnan(target_value)
                error('Invalid input. Please enter a valid expression.');
            end

            % Define the target value as a function of time if it contains "t"
            if contains(target_value_str, 't')
                target_value = str2func(['@(t)', target_value_str]);
            else
                % If it's a constant value, create a function that returns the constant + 0*t
                target_value = str2func(['@(t)', target_value_str, ' + 0*t']);
            end
        catch
            uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
            target_value = [];
        end
    end

    % Display the target value type
    if ~isempty(target_value)
        disp('The target value has been set.');
    end
end
