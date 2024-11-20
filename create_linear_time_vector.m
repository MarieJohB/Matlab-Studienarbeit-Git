function [t] = create_linear_time_vector(num_sections)
%create a linear time vector t
% with the input elementes: vector with start times, vector with end times
% and vector with time steps 

% call function to get time sections
[start_time, end_time, time_steps] = get_time_vector(num_sections);



min_start_time = min(start_time);
max_end_time = max(end_time);

% check whether the target value was defined piecewise or continuous
num_sections = length(end_time);


if num_sections == 1 % target value is continuous
        % Define the time vector
        num_steps = round((max_end_time - min_start_time) / min(time_steps)); % calculate integer number of steps
        t = linspace(start_time, end_time, num_steps + 1)'; % create time vector with specified values, transpose to column vector


    elseif num_sections > 1 % target value is defined piecewise
        % Initialize the time vector and target array
        t = [];
        
        
        for i = 1:num_sections
            % useing min(time_steps) for all because when using lsim: time has to
            % be spaced venly
            num_steps = round((end_time(i) - start_time(i)) / min(time_steps)); % calculate integer number of steps
            t_section = linspace(start_time(i), end_time(i), num_steps + 1)'; % create equal intervals for each section with time_steps(i)
            
            % Avoid duplicate values at the transition points
            % neccessary because end_time(1) = start_time(2) and so on 
            if i > 1
                t_section = t_section(2:end); % remove the first element to avoid duplicate
            end

            t = [t; t_section]; % concatenate as column vector
            % and transpond t so its applicable for lsim
        end
end

end