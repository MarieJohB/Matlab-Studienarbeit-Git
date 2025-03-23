function [t] = create_linear_time_vector_disturbance(num_sections, disturbance_label)
% CREATE_LINEAR_TIME_VECTOR_DISTURBANCE - Creates a linear time vector for a disturbance signal
% Similar to create_linear_time_vector but considers the specific time parameters of the disturbance
%
% Parameters:
%   num_sections - Number of time sections
%   disturbance_label - Label for the disturbance (e.g., 'd_1' or 'd_2')
%
% Returns:
%   t - Combined linear time vector for all sections

% Initialize time vector
t = [];

% Check if num_sections is valid
if isempty(num_sections)
    disp('Operation cancelled: No sections defined.');
    return;
end

% Get time parameters
[start_time, end_time, time_steps] = get_time_vector_disturbance(num_sections, [], [], [], disturbance_label);

% Check if time parameters were returned
if isempty(start_time) || isempty(end_time) || isempty(time_steps)
    disp('Operation cancelled: Time parameters not provided.');
    return;
end

% Create time vector for each section
for i = 1:num_sections
    % Calculate number of points in this section
    num_points = round((end_time(i) - start_time(i)) / time_steps(i)) + 1;
    
    % Create time vector for this section
    section_t = linspace(start_time(i), end_time(i), num_points);
    
    % If this is not the first section, remove the first point to avoid duplication
    if i > 1 && ~isempty(t) && ~isempty(section_t)
        section_t = section_t(2:end);
    end
    
    % Append to overall time vector
    t = [t, section_t];
end

% Display time vector information
if ~isempty(t)
    disp(['Created ', disturbance_label, ' time vector with ', num2str(length(t)), ' points from ', ...
        num2str(t(1)), ' to ', num2str(t(end)), ' across ', num2str(num_sections), ' section(s).']);
else
    disp(['Warning: Created an empty time vector for ', disturbance_label, '.']);
end
end