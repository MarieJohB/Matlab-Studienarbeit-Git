function disturbance = get_disturbance_input_version2(start_time, end_time, time_steps)


num_sections = length(start_time);

min_start_time = min(start_time);
max_end_time = max(end_time);
min_time_step = min(time_steps);

    if num_sections == 1
    % 1 section: function r(t) is continous
    % also possible/not yet developed: disturbance is defined piecewise although r(t) is
    % continous

    message = sprintf('The disturbance will be defined from start time= %.2f until end time %.2f \n Please enter an expression for the disturbance', min_start_time, max_end_time);
    uiwait(msgbox(message, 'Information', 'modal'));

    num_steps = (end_time - start_time) / min_time_step;
    t = linspace(min_start_time, max_end_time, num_steps + 1);

    disturbance = get_function_input();


    elseif num_sections > 1
    % more than one section: function r(t) is defined piecewise
    % user can decide whether the sections shall be kept or disturbance shall
    % be continous as well
    % also possible/not yet developed: disturbance has different sections than
    % r(t)


    % ask user 
    user_decision = button_dialog_disturbance();


    if user_decision == true
        % disturbances are defined piecewise
        for i = 1:num_sections 
            disturbance{i} = get_function_input;
        end

    elseif user_decision == false
        % disturbance is defined continous
        message = sprintf('The disturbance will be defined from start time= %.2f until end time %.2f \n Please enter an expression for the disturbance', min_start_time, max_end_time);
        uiwait(msgbox(message, 'Information', 'modal'));
        num_steps = (end_time - start_time) / min_time_step;
        t = linspace(min_start_time, max_end_time, num_steps + 1);

        disturbance = get_function_input();
    end


    elseif num_scetions < 1 
        disp('Error: number of sections is smaller than 1!');
    
    
    end

end




function user_decision = button_dialog_disturbance()

   %create figure to ask user
    fig = uifigure('Position', [100, 100, 300, 150], 'Name', 'Disturbance continuous or piecewise');

     uilabel(fig, 'Position', [5, 110, 300, 50], 'Text', ... 
         'Is the disturbance continous or piecewise defined:', ... 
         'FontSize', 12, 'HorizontalAlignment', 'center');
    
    % initialze variable for flag
    user_choice = false;

        btn1 = uibutton(fig, 'Position', [50, 50, 100, 50], 'Text', 'Continuous', ...
        'ButtonPushedFcn', @(btn1,event) set_and_close(false));

        btn2 = uibutton(fig, 'Position', [150, 50, 100, 50], 'Text', 'Piecewise', ...
        'ButtonPushedFcn', @(btn2,event) set_and_close(true));

    uiwait(fig);
    
    % set the user choise
    function set_and_close(choice)
        user_choice = choice;
        uiresume(fig);
        close(fig); % close figure
    end
user_decision = user_choice;
end
