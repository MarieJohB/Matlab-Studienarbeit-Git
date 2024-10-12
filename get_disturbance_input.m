function disturbance = get_disturbance_input()
    % Initialize disturbance to be empty
    disturbance = [];
    
    % Loop until valid input is provided or user cancels the dialog
    while isempty(disturbance)
        % Prompt for disturbance value (constant or time-dependent)
        prompt = {'Enter the disturbance value (use "t" for time-dependent values, e.g., "5" or "sin(t)"): '};
        dlgtitle = 'Disturbance Value Input';
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
            disturbance_str = strrep(answer{1}, ',', '.');

            % Test the input
            t = 0; % Dummy value for validation
            disturbance = eval(disturbance_str);
            if isnan(disturbance)
                error('Invalid input. Please enter a valid expression.');
            end

            % Define the disturbance value as a function of time if it contains "t"
            if contains(disturbance_str, 't')
                disturbance = str2func(['@(t)', disturbance_str]);
            else
                % If it's a constant value, create a function that returns the constant + 0*t
                disturbance = str2func(['@(t)', disturbance_str, ' + 0*t']);
            end
        catch
            uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
            disturbance = [];
        end
    end

    % Display the disturbance type
    if ~isempty(disturbance)
        disp('The disturbance value has been set.');
    end
end