function [start_time, end_time, time_steps] = get_time_vector(num_sections)

% initializing persistent variables to check if the function was called
% before and the time vector was already entered be the user
persistent start_value;
persistent end_value;
persistent steps_value;

if isempty(start_value) || isempty(end_value) ||isempty(steps_value)
% if the persistend variables are empty: fill them 
 % initalize variables to for  start time, end time and steps
 start_time = NaN * ones(1, num_sections);
 end_time = NaN * ones(1, num_sections);
 time_steps = NaN * ones(1, num_sections);
        


 for i = 1:num_sections
    
        while isnan(start_time(i)) || isnan(end_time(i)) || isnan(time_steps(i))
            
            if i==1
             
            % define request for user input
            prompt = {'Start Time:', 'End Time:', 'Time Steps:'};
            dlgtitle = 'Input';
            dims = [1 35];
            answer = inputdlg(prompt, dlgtitle, dims);

                % Check if user cancelled the dialog
                if isempty(answer)
                    disp('Operation cancelled by user.');
                       return;
                end

            start_time_str = strrep(answer{1}, ',', '.');
            start_time(i) = str2double(start_time_str);

            end_time_str = strrep(answer{2}, ',', '.');
            end_time(i) = str2double(end_time_str);

            time_steps_str = strrep(answer{3}, ',', '.');
            time_steps(i) = str2double(time_steps_str);

            % check if end and start > 0;

            %checking if end > start:
            if start_time(i) >= end_time(i)
                uiwait(msgbox('Invalid input. Please make sure that end happens after start.', 'Error', 'error'));
                start_time(i) = NaN;
                end_time(i) = NaN;
            end

            % check if end and start > 0;
            if start_time(i) < 0
                uiwait(msgbox('Invalid input. Please enter positive values only.', 'Error', 'error'));
                start_time(i) = NaN;
            end

            % checking if step size is compatible with time span
            numSteps = (end_time(i) - start_time(i)) / time_steps(i);
            if mod(numSteps, 1) ~= 0
                uiwait(msgbox('Invalid input. Please enter suitable size for time steps.', 'Error', 'error'));
                time_steps(i) = NaN;
            end

            
            elseif i > 1 
                % for the following sections
                % end time is the start time for the next section
                start_time(i) = end_time(i-1);

                %adding information for user: start time of current section
                message = sprintf('Start time of this section = %f. \nPlease enter end time and steps for this section', start_time(i));
                uiwait(msgbox(message, 'Information', 'modal'));
                % define request for user input
                prompt = {'End Time:', 'Time Steps:'};
                dlgtitle = 'Input';
                dims = [1 35];
 
                answer = inputdlg(prompt, dlgtitle, dims);

               
                % Check if user cancelled the dialog
                if isempty(answer)
                    disp('Operation cancelled by user.');
                       return;
                end



                end_time_str = strrep(answer{1}, ',', '.');
                end_time(i) = str2double(end_time_str);
    
                time_steps_str = strrep(answer{2}, ',', '.');
                time_steps(i) = str2double(time_steps_str);

                %checking if end > start:
                if start_time(i) > end_time(i)
                    uiwait(msgbox('Invalid input. Please make sure that end happens after start.', 'Error', 'error'));
                    start_time(i) = NaN;
                    end_time(i) = NaN;
                end

                % checking if step size is compatible with time span
                numSteps = (end_time(i) - start_time(i)) / time_steps(i);
                if mod(numSteps, 1) ~= 0
                    uiwait(msgbox('Invalid input. Please enter suitable size for time steps.', 'Error', 'error'));
                    time_steps(i) = NaN;
                end

             end
         end

% save the defined values in persistend variables
start_value = start_time;
end_value = end_time;
steps_value = time_steps;

 end




else % if time vector has been set already
% saved values in the persistent variables are given back as the previously
% defined time vector

start_time = start_value;
end_time = end_value;
time_steps = steps_value;


end
end