%---------------------------------------------------------------------
function [y_inf_1, e_inf_1, y_inf_2, e_inf_2, e_inf_3, cansysjump] = calculateStationaryValues(T_s, S_s, R_s, D1_s, D2_s, G_s)
    % calculateStationaryValues calculates the steady-state values of the system
    % via simulation. It supports continuous as well as piecewise inputs for R, D1, and D2.
    %
    % If T_s, S_s, or G_s are numeric, they are converted to tf objects (ensuring row vectors).
    % The steady-state values are obtained by simulating the step response over a fixed time interval
    % (t_final) and taking the final value.
    %
    % Additionally, the "jumpability" is estimated by evaluating T(s)*R(s) at a high frequency.
    
    % Ensure T_s, S_s, and G_s are tf objects
    if ~isa(T_s, 'tf')
        T_s = tf(T_s(:)', 1);
    end
    if ~isa(S_s, 'tf')
        S_s = tf(S_s(:)', 1);
    end
    if ~isa(G_s, 'tf')
        G_s = tf(G_s(:)', 1);
    end
    
    % Use a default simulation horizon (adjust as needed)
    t_final = 100;
    
    % Helper: Ensure input is a cell array
    function cellOut = ensureCell(x)
        if iscell(x)
            cellOut = x;
        else
            cellOut = {x};
        end
    end

    % Wrap R_s, D1_s, and D2_s in cell arrays if needed
    R_cell  = ensureCell(R_s);
    D1_cell = ensureCell(D1_s);
    D2_cell = ensureCell(D2_s);
    
    % Determine the number of segments (nSeg) from the maximum length of the cell arrays
    nSeg = max([numel(R_cell), numel(D1_cell), numel(D2_cell)]);
    
    % Preallocate outputs (one value per segment)
    y_inf_1   = zeros(1, nSeg);
    e_inf_1   = zeros(1, nSeg);
    y_inf_2   = zeros(1, nSeg);
    e_inf_2   = zeros(1, nSeg);
    e_inf_3   = zeros(1, nSeg);
    cansysjump = zeros(1, nSeg);
    
    % Helper: Select an element from a cell array; if missing, use the first element.
    function out = selectSegment(cellArray, idx)
        if numel(cellArray) >= idx
            out = cellArray{idx};
        else
            out = cellArray{1};
        end
    end

    % Loop over each segment
    for i = 1:nSeg
        % Select segments for R, D1, D2
        R_i  = selectSegment(R_cell, i);
        D1_i = selectSegment(D1_cell, i);
        D2_i = selectSegment(D2_cell, i);
        
        % Convert R_i, D1_i, and D2_i to tf objects if necessary (and ensure row vectors)
        if ~isa(R_i, 'tf')
            if isnumeric(R_i)
                R_i = tf(R_i(:)', 1);
            else
                error('R_i must be numeric or a tf object.');
            end
        end
        if ~isa(D1_i, 'tf')
            if isnumeric(D1_i)
                D1_i = tf(D1_i(:)', 1);
            else
                error('D1_i must be numeric or a tf object.');
            end
        end
        if ~isa(D2_i, 'tf')
            if isnumeric(D2_i)
                D2_i = tf(D2_i(:)', 1);
            else
                error('D2_i must be numeric or a tf object.');
            end
        end
        
        % Compute steady-state values via simulation (using step response)
        % y_inf_1: final output of T_s * R_i
        sys_y1 = T_s * R_i;
        [y_sim, ~] = step(sys_y1, t_final);
        y_inf_1(i) = y_sim(end);
        
        % e_inf_1: final value of S_s * R_i
        sys_e1 = S_s * R_i;
        [y_sim, ~] = step(sys_e1, t_final);
        e_inf_1(i) = y_sim(end);
        
        % y_inf_2: final value of S_s * D1_i
        sys_y2 = S_s * D1_i;
        [y_sim, ~] = step(sys_y2, t_final);
        y_inf_2(i) = y_sim(end);
        
        % e_inf_2: final value of -S_s * D1_i
        sys_e2 = -S_s * D1_i;
        [y_sim, ~] = step(sys_e2, t_final);
        e_inf_2(i) = y_sim(end);
        
        % e_inf_3: final value of S_s * R_i - G_s * S_s * D2_i
        sys_e3 = S_s * R_i - G_s * S_s * D2_i;
        [y_sim, ~] = step(sys_e3, t_final);
        e_inf_3(i) = y_sim(end);
        
        % Calculate jumpability as the high-frequency gain of T_s*R_i
        try
            hf = 1e6;
            TR_val = evalfr(sys_y1, 1i*hf);
            cansysjump(i) = abs(TR_val * hf);
        catch
            cansysjump(i) = NaN;
        end
    end
    
    % If only one segment exists, return scalar values
    if nSeg == 1
        y_inf_1   = y_inf_1(1);
        e_inf_1   = e_inf_1(1);
        y_inf_2   = y_inf_2(1);
        e_inf_2   = e_inf_2(1);
        e_inf_3   = e_inf_3(1);
        cansysjump = cansysjump(1);
    end
end
