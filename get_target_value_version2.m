
function target_value = get_target_value_version2()
    % output: 
    % if continous: 
    % Initialize target_value to be empty
    target_value = [];

    % Create a figure for user interaction
    selection_fig = uifigure('Name', 'Select Function Type', 'Position', [500, 500, 300, 150]);
    
    % Add buttons for continuous and piecewise functions
    uibutton(selection_fig, 'Position', [50, 80, 200, 50], 'Text', 'Continuous Function', ...
        'ButtonPushedFcn', @(btn, event) continuous_function_input(selection_fig));
    uibutton(selection_fig, 'Position', [50, 20, 200, 50], 'Text', 'Define Function in Sections', ...
        'ButtonPushedFcn', @(btn, event) piecewise_function_input(selection_fig));

    % continue in the main function after target_value was entered in other
    % functions
    uiwait(selection_fig);

    % set the target_value in the main function 
    target_value = getappdata(selection_fig, 'target_value');

    % Close the selection figure
    close(selection_fig);
end


function continuous_function_input(selection_fig)

    
    % get the users input
    target_value = get_function_input();

    % give the target_function to the main function
    setappdata(selection_fig, 'target_value', target_value);

     % Resume execution in the main function 
     uiresume(selection_fig);

end

function piecewise_function_input(selection_fig)
    
    
    % initalize variable for number of sections
    num_sections = [];

    while isempty(num_sections)
    prompt = {'Enter the number of sections (max 5): '};
    dlgtitle = 'Number of Sections';
    dims = [1 50];
    answer = inputdlg(prompt, dlgtitle, dims);
    
    % Check if user cancelled the dialog
    if isempty(answer)
        disp('Operation cancelled by user.');
        return;
    end
    
    % Get the number of sections
    % maximum number of sections is 5
    % number cannot be decimal

    % Replace commas with periods
    num_sections_str = strrep(answer{1}, ',', '.');
    num_sections = str2double(num_sections_str);

    if isnan(num_sections) || num_sections < 1 || num_sections > 5 || mod(num_sections, 1) ~= 0
        uiwait(msgbox('Invalid input. Please enter a number between 1 and 5.', 'Error', 'error'));
        num_sections = [];
    end

    end

    % get time vector
    [start_time, end_time, time_steps] = get_time_vector(num_sections);

    % get function for every time section 
    % initialize array for functions
    target_value = cell(1, num_sections);
% fill array with functions created by user
for i = 1:num_sections

    %adding information for user: time sections
    message = sprintf('The following time vector has been set for section %.0f: \nStart Time: %.2f \nEnd Time: %.2f \nTime Steps: %.2f \nPlease enter a function', i, start_time(i), end_time(i), time_steps(i));
    uiwait(msgbox(message, 'Information', 'modal'));

     target_value{i} = get_function_input;
end

% show function
disp('Stored functions:');
disp(target_value);


disp('Start times:'); 
disp(start_time);
disp('End times:'); 
disp(end_time);
disp('Time steps:'); 
disp(time_steps);
disp('Target values:'); 
disp(target_value);

 % give the target_function to the main function
  setappdata(selection_fig, 'target_value', target_value);
 % Resume execution in the main function 
  uiresume(selection_fig);

end





