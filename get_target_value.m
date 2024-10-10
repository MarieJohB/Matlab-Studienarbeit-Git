function target_value = get_target_value()
    % Prompt user to choose target type
    target_type = menu('Choose the type of target value:', 'Constant', 'Time-Dependent');

    switch target_type
        case 1 % Constant target value
            prompt = {'Enter the constant target value:'};
            dlgtitle = 'Constant Target Value Input';
            dims = [1 50];
            answer = inputdlg(prompt, dlgtitle, dims);
            
            % Validate input
            try
                target_value = str2double(answer{1});
                if isnan(target_value)
                    error('Invalid input. Please enter a numeric value.');
                end
            catch
                uiwait(msgbox('Invalid input. Please enter the value correctly.', 'Error', 'error'));
                target_value = get_target_value(); % Recursively call the function until a valid input is received
            end

        case 2 % Time-dependent target value
            prompt = {'Enter the expression for the time-dependent target value (use "t" as the variable):'};
            dlgtitle = 'Time-Dependent Target Value Input';
            dims = [1 50];
            answer = inputdlg(prompt, dlgtitle, dims);
            
            % Validate input
            try
                % Test the input with a dummy value of t
                t = 0;
                target_value = eval(answer{1});
                if isnan(target_value)
                    error('Invalid input. Please enter a valid expression.');
                end
                % Define the target value as a function of time
                target_value = str2func(['@(t)' answer{1}]);
            catch
                uiwait(msgbox('Invalid input. Please enter the expression correctly.', 'Error', 'error'));
                target_value = get_target_value(); % Recursively call the function until a valid input is received
            end

        otherwise
            disp('No target type selected.');
            target_value = [];
            return;
    end

    % Display the target value type
    if ~isempty(target_value)
        if target_type == 1
            disp(['The constant target value is: ', num2str(target_value)]);
        else
            disp('The time-dependent target value has been set.');
        end
    end
end