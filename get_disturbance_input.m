function disturbance = get_disturbance_input(name)
    % Prompt user to choose disturbance type
    disturbance_type = menu(['Choose the type of disturbance ', name, ':'], 'Constant', 'Time-Dependent');

    switch disturbance_type
        case 1 % Constant disturbance
            prompt = {['Enter the constant disturbance value for ', name, ':']};
            dlgtitle = [name, ' Constant Disturbance Input'];
            dims = [1 50];
            answer = inputdlg(prompt, dlgtitle, dims);

            % Validate input
            try
                disturbance = str2double(answer{1});
                if isnan(disturbance)
                    error('Invalid input. Please enter a numeric value.');
                end
                disturbance = disturbance * ones(1, 1000); % Convert to time vector
            catch
                uiwait(msgbox('Invalid input. Please enter the value correctly.', 'Error', 'error'));
                disturbance = get_disturbance_input(name); % Recursively call the function until a valid input is received
            end

        case 2 % Time-dependent disturbance
            prompt = {['Enter the expression for the time-dependent disturbance ', name, ' (use "t" as the variable):']};
            dlgtitle = [name, ' Time-Dependent Disturbance Input'];
            dims = [1 50];
            answer = inputdlg(prompt, dlgtitle, dims);

            % Validate input
            try
                % Test the input with a dummy value of t
                t = 0;
                disturbance_value = eval(answer{1});
                if isnan(disturbance_value)
                    error('Invalid input. Please enter a valid expression.');
                end
                % Define the disturbance value as a function of time
                disturbance_func = str2func(['@(t)', answer{1}]);
                t = linspace(0, 10, 1000); % Time vector
                disturbance = arrayfun(disturbance_func, t); % Evaluate the function handle over time
            catch
                uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
                disturbance = get_disturbance_input(name); % Recursively call the function until a valid input is received
            end

        otherwise
            disp('No valid disturbance type selected.');
            disturbance = [];
            return;
    end

    % Display the disturbance type
    if ~isempty(disturbance)
        if disturbance_type == 1
            disp(['The constant disturbance value for ', name, ' is: ', num2str(disturbance(1))]);
        else
            disp(['The time-dependent disturbance ', name, ' has been set.']);
        end
    end
end