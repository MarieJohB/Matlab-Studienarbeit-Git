function [d1, d2, t] = get_disturbances()

% function is not used in app?! 
    % Define the time vector
    t = linspace(0, 10, 1000); % 0 to 10 seconds with 1000 points

    % Get D1 disturbance
    d1 = get_disturbance_input();

    % Get D2 disturbance
    d2 = get_disturbance_input();
    
    
end