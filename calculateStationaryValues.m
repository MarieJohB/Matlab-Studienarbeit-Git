function [y_stat_TR, e_stat_SR, y_stat_SD1, e_stat_SD1, e_stat, cansysjump, inputs_constant] = calculateStationaryValues(T, S, r_s, d_1_s, d_2_s, G, r, d_1, d_2, t, input_r, input_d1)
    % CALCULATESTATIONARYVALUES - Calculate steady-state values for a control system
    % This function calculates steady-state values primarily using simulation results
    %
    % Inputs:
    %   T       - Complementary sensitivity function T(s) = L(s)/(1+L(s))
    %   S       - Sensitivity function S(s) = 1/(1+L(s))
    %   r_s     - Reference input value (scalar for step input)
    %   d_1_s   - Output disturbance value (scalar for step input)
    %   d_2_s   - Input disturbance value (scalar for step input)
    %   G       - Plant transfer function G(s)
    %   r       - Reference input function(s)
    %   d_1     - Output disturbance function(s)
    %   d_2     - Input disturbance function(s)
    %   t       - Time vector for simulation
    %   input_r - Reference input values corresponding to time vector t
    %   input_d1- Disturbance input values corresponding to time vector t
    %
    % Outputs:
    %   y_stat_TR     - Steady-state output for reference tracking only
    %   e_stat_SR     - Steady-state error for reference tracking only
    %   y_stat_SD1    - Steady-state output for output disturbance only
    %   e_stat_SD1    - Steady-state error for output disturbance only
    %   e_stat        - Total steady-state error with all inputs
    %   cansysjump    - Flag indicating if system is "jumpable" (1 if yes, 0 if no)
    %   inputs_constant - Structure with flags indicating if inputs are constant
    
    % Initialize outputs
    y_stat_TR = NaN;
    e_stat_SR = NaN;
    y_stat_SD1 = NaN;
    e_stat_SD1 = NaN;
    e_stat = NaN;
    cansysjump = 0;
    
    % Initialize inputs_constant structure
    inputs_constant = struct('reference', true, 'disturbance1', true, 'disturbance2', true);
    
    % Check if inputs are constant
    if nargin >= 7
        inputs_constant.reference = isInputConstant(r);
    end
    if nargin >= 8
        inputs_constant.disturbance1 = isInputConstant(d_1);
    end
    if nargin >= 9
        inputs_constant.disturbance2 = isInputConstant(d_2);
    end
    
    % Simulation-based calculation (primary method)
    if nargin >= 11 && ~isempty(t) && ~isempty(input_r)
        % Setpoint response via simulation
        try
            y_TR = lsim(T, input_r, t);
            e_SR = lsim(S, input_r, t);
            
            if ~isempty(y_TR) && length(y_TR) > 10  % Ensure enough points for reliable estimate
                % Use the last 10% of simulation points for more reliable steady-state value
                endpart = max(1, round(0.9*length(y_TR)));
                y_stat_TR = mean(y_TR(endpart:end));
                e_stat_SR = mean(e_SR(endpart:end));
                
                disp(['Simulation-based calculation - Output: ', num2str(y_stat_TR), ', Error: ', num2str(e_stat_SR)]);
            end
        catch ME
            disp(['Simulation error for reference tracking: ', ME.message]);
        end
    end
    
    % Disturbance response via simulation
    if nargin >= 12 && ~isempty(t) && ~isempty(input_d1)
        try
            y_SD1 = lsim(-S, input_d1, t);
            e_SD1 = lsim(S, input_d1, t);
            
            if ~isempty(y_SD1) && length(y_SD1) > 10
                % Use the last 10% of simulation points for more reliable steady-state value
                endpart = max(1, round(0.9*length(y_SD1)));
                y_stat_SD1 = mean(y_SD1(endpart:end));
                e_stat_SD1 = mean(e_SD1(endpart:end));
                
                disp(['Simulation-based calculation - Disturbance: ', num2str(y_stat_SD1), ', Error: ', num2str(e_stat_SD1)]);
            end
        catch ME
            disp(['Simulation error for disturbance response: ', ME.message]);
        end
    end
    
    % If simulation didn't provide values, try analytical methods as fallback
    if isnan(y_stat_TR) || isnan(e_stat_SR)
        disp('Simulation-based calculation failed. Using analytical fallback...');
        
        % For reference tracking only (analytical fallback)
        if isscalar(r_s) && inputs_constant.reference && ~isempty(T) && ~isnan(T)
            try
                % Method 1: DC gain calculation
                [num_T, den_T] = tfdata(T, 'v');
                if den_T(end) ~= 0  % No pole at s=0
                    % Direct DC gain calculation
                    dc_gain_T = num_T(end) / den_T(end);
                    y_stat_TR = dc_gain_T * r_s;
                    
                    [num_S, den_S] = tfdata(S, 'v');
                    if den_S(end) ~= 0  % No pole at s=0
                        dc_gain_S = num_S(end) / den_S(end);
                        e_stat_SR = dc_gain_S * r_s;
                    else
                        % Integrator in the loop
                        e_stat_SR = 0;
                    end
                else
                    % Method 2: Using small s-value
                    s_small = 1e-9;
                    y_stat_TR = evalfr(T, s_small) * r_s;
                    e_stat_SR = evalfr(S, s_small) * r_s;
                end
                
                disp(['Analytical calculation - Output: ', num2str(y_stat_TR), ', Error: ', num2str(e_stat_SR)]);
            catch ME
                disp(['Analytical calculation error: ', ME.message]);
            end
        end
        
        % For output disturbance only (analytical fallback)
        if isnan(y_stat_SD1) || isnan(e_stat_SD1)
            if isscalar(d_1_s) && inputs_constant.disturbance1 && ~isempty(S) && ~isnan(S)
                try
                    % Method 1: DC gain calculation
                    [num_S, den_S] = tfdata(S, 'v');
                    if den_S(end) ~= 0  % No pole at s=0
                        dc_gain_S = num_S(end) / den_S(end);
                        y_stat_SD1 = -dc_gain_S * d_1_s;
                        e_stat_SD1 = dc_gain_S * d_1_s;
                    else
                        % Method 2: Using small s-value
                        s_small = 1e-9;
                        y_stat_SD1 = -evalfr(S, s_small) * d_1_s;
                        e_stat_SD1 = evalfr(S, s_small) * d_1_s;
                    end
                    
                    disp(['Analytical calculation - Disturbance: ', num2str(y_stat_SD1), ', Error: ', num2str(e_stat_SD1)]);
                catch ME
                    disp(['Analytical calculation error: ', ME.message]);
                end
            end
        end
    end
    
    % Calculate total steady-state error
    if ~isnan(e_stat_SR) && ~isnan(e_stat_SD1)
        e_stat = e_stat_SR + e_stat_SD1;
    elseif ~isnan(e_stat_SR)
        e_stat = e_stat_SR;
    elseif ~isnan(e_stat_SD1)
        e_stat = e_stat_SD1;
    end
    
    % Check if system is "jumpable" (zero steady-state error for step input)
    % First check based on simulation results
    if ~isnan(e_stat_SR) && abs(e_stat_SR) < 1e-3
        cansysjump = 1;  % System is jumpable based on small error
    else
        % Check for integrator in the loop
        try
            % Extract open-loop transfer function
            [num_S, den_S] = tfdata(S, 'v');
            if any(abs(den_S) < 1e-10)
                cansysjump = 1;  % System has pole at origin (integrator)
            else
                % Use a frequency response approach
                s_small = 1e-9;
                mag_S_0 = abs(evalfr(S, s_small));
                if mag_S_0 < 1e-3
                    cansysjump = 1;  % Very small sensitivity at DC indicates integrator
                else
                    cansysjump = 0;  % No integrator detected
                end
            end
        catch
            % Default to checking simulation error
            cansysjump = (~isnan(e_stat_SR) && abs(e_stat_SR) < 1e-3);
        end
    end
end

function isConstant = isInputConstant(input)
    % Check if the input is constant (scalar, empty, or contains constant function)
    isConstant = true;
    
    % If it's a scalar, it's constant
    if isscalar(input)
        return;
    end
    
    % If it's empty, consider it constant
    if isempty(input)
        return;
    end
    
    % If it's a cell array (piecewise function), check each piece
    if iscell(input)
        for i = 1:length(input)
            if isa(input{i}, 'function_handle')
                % Check if function contains 't' in its expression
                func_str = func2str(input{i});
                if contains(func_str, 't') && ~contains(func_str, '0*t')
                    isConstant = false;
                    return;
                end
            end
        end
    elseif isa(input, 'function_handle')
        % For a single function handle, check if it contains 't'
        func_str = func2str(input);
        if contains(func_str, 't') && ~contains(func_str, '0*t')
            isConstant = false;
            return;
        end
    end
end