function [y_inf_1, e_inf_1, y_inf_2, e_inf_2, e_inf_3, cansysjump] = calculateStationaryValues(T_s, S_s, R_s, D1_s, D2_s, G_s)
    % calculateStationaryValues calculates the stationary values of the control loop
    % the function supports continous functions as well as piecewise
    % defined functions for r, d1 and d2
    % should the inputs be of format double, a conversion into tf-objects
    % is performed 
    % It is also made sure that numeric inputs are defined as row-vectors
    % to avoid issues with tf
    %
    % \ac{Robustheit} und \ac{Flexibilität} werden durch die Unterstützung
    % unterschiedlicher Eingabetypen erreicht
    % \cite{MatlabBestPractices2019}.???
    
    % ensuring T_s, S_s and G_s are tf-objects
    if ~isa(T_s, 'tf')
        T_s = tf(T_s, 1);
    end
    if ~isa(S_s, 'tf')
        S_s = tf(S_s, 1);
    end
    if ~isa(G_s, 'tf')
        G_s = tf(G_s, 1);
    end
    
    % convert T_s, S_s and G_s to symbolic transfer functions
    % use of s = sym('s') istead of syms s, to avoid problems with static
    % work spaces
    s = sym('s');
    [T_num, T_den] = tfdata(T_s, 'v');
    [S_num, S_den] = tfdata(S_s, 'v');
    [G_num, G_den] = tfdata(G_s, 'v');
    
    T_sys = poly2sym(T_num, s) / poly2sym(T_den, s);
    S_sys = poly2sym(S_num, s) / poly2sym(S_den, s);
    G_sys = poly2sym(G_num, s) / poly2sym(G_den, s);
    
    % help-function: ensuring input to be a cell array
    function cellOut = ensureCell(x)
        if iscell(x)
            cellOut = x;
        else
            cellOut = {x};
        end
    end

    % check R_s, D1_s and D2_s are cell arrays
    % import to support continous functions and piecewise
    % defined functions
    R_cell  = ensureCell(R_s);
    D1_cell = ensureCell(D1_s);
    D2_cell = ensureCell(D2_s);
    
    % find number of segments using longest length cell array
    nSeg = max([numel(R_cell), numel(D1_cell), numel(D2_cell)]);
    
    % Initialization
    y_inf_1   = zeros(1, nSeg);
    e_inf_1   = zeros(1, nSeg);
    y_inf_2   = zeros(1, nSeg);
    e_inf_2   = zeros(1, nSeg);
    e_inf_3   = zeros(1, nSeg);
    cansysjump = zeros(1, nSeg);
    
    %help-function for selection of an element: for the case that segment i
    %does not exsist, the firt element is used (continous input)
    function out = selectSegment(cellArray, idx)
        if numel(cellArray) >= idx
            out = cellArray{idx};
        else
            out = cellArray{1};
        end
    end

    % loop over the segments 
    for i = 1:nSeg
        % selecting the segments (oder continous values, in case of no array)
        R_i  = selectSegment(R_cell, i);
        D1_i = selectSegment(D1_cell, i);
        D2_i = selectSegment(D2_cell, i);
        
        % conversion of inputs to tf-objects if neccessary. 
        % Ensuring numeric inputs are row vectors
        if ~isa(R_i, 'tf')
            if isnumeric(R_i)
                R_i = tf(R_i(:)', 1);
            else
                error('R_i must be numeric oder or a tf-object.');
            end
        end
        if ~isa(D1_i, 'tf')
            if isnumeric(D1_i)
                D1_i = tf(D1_i(:)', 1);
            else
                error('D1_i must be numeric oder or a tf-object.');
            end
        end
        if ~isa(D2_i, 'tf')
            if isnumeric(D2_i)
                D2_i = tf(D2_i(:)', 1);
            else
                error('D2_i must be numeric oder or a tf-object.');
            end
        end
        
        % Conversion onto symbolic transfer function
        [R_num, R_den]   = tfdata(R_i, 'v');
        [D1_num, D1_den] = tfdata(D1_i, 'v');
        [D2_num, D2_den] = tfdata(D2_i, 'v');
        
        R_sys  = poly2sym(R_num, s) / poly2sym(R_den, s);
        D1_sys = poly2sym(D1_num, s) / poly2sym(D1_den, s);
        D2_sys = poly2sym(D2_num, s) / poly2sym(D2_den, s);
        
        % calculating the stationary values for the current segment
        y_inf_1(i)    = double(limit(s * T_sys * R_sys, s, 0));
        e_inf_1(i)    = double(limit(s * S_sys * R_sys, s, 0));
        y_inf_2(i)    = double(limit(s * S_sys * D1_sys, s, 0));
        e_inf_2(i)    = double(limit(-s * S_sys * D1_sys, s, 0));
        e_inf_3(i)    = double(limit(s * (S_sys * R_sys - G_sys * S_sys * D2_sys), s, 0));
        cansysjump(i) = double(limit(s * T_sys * R_sys, s, inf));
    end

    % For the case of only 1 segment, scalar values are returned
    if nSeg == 1
        y_inf_1    = y_inf_1(1);
        e_inf_1    = e_inf_1(1);
        y_inf_2    = y_inf_2(1);
        e_inf_2    = e_inf_2(1);
        e_inf_3    = e_inf_3(1);
        cansysjump = cansysjump(1);
    end
end
