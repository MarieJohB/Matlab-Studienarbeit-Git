function [num_c, den_c] = conditionTransferFunction(num, den)
    % CONDITIONTRANSFERFUNCTION Improve numerical conditioning of transfer function
    % Applies scaling and cleaning to improve numerical properties
    %
    % Inputs:
    %   num - Numerator coefficients vector
    %   den - Denominator coefficients vector
    %
    % Outputs:
    %   num_c - Conditioned numerator coefficients
    %   den_c - Conditioned denominator coefficients
    
    % 1. Remove leading zeros
    num = removeLeadingZeros(num);
    den = removeLeadingZeros(den);
    
    % 2. Remove very small coefficients (numerical noise)
    num_c = cleanSmallCoefficients(num);
    den_c = cleanSmallCoefficients(den);
    
    % 3. Apply scaling for better numerical conditioning
    [num_c, den_c] = applyScaling(num_c, den_c);
    
    % 4. Round very small numbers to zero
    num_c = roundToZero(num_c, 1e-12);
    den_c = roundToZero(den_c, 1e-12);
    
    % 5. Normalize to make leading coefficient of denominator 1
    [num_c, den_c] = normalizeCoefficients(num_c, den_c);
end

function coeff = removeLeadingZeros(coeff)
    % Remove leading zeros from coefficient vector
    idx = find(abs(coeff) > eps, 1, 'first');
    
    if isempty(idx)
        coeff = 0;  % All zeros
    else
        coeff = coeff(idx:end);
    end
end

function coeff = cleanSmallCoefficients(coeff)
    % Remove small coefficients that are likely numerical noise
    
    % Determine threshold based on coefficient magnitude
    max_coeff = max(abs(coeff));
    threshold = max_coeff * 1e-10;
    
    % Round small values to zero
    coeff(abs(coeff) < threshold) = 0;
end

function [num, den] = applyScaling(num, den)
    % Apply scaling to improve numerical conditioning
    
    % Get maximum coefficient values
    max_num = max(abs(num));
    max_den = max(abs(den));
    
    % Get minimum non-zero coefficient values
    min_num = min(abs(num(abs(num) > eps)));
    min_den = min(abs(den(abs(den) > eps)));
    
    if isempty(min_num), min_num = max_num; end
    if isempty(min_den), min_den = max_den; end
    
    % Check dynamic range of coefficients
    num_range = max_num / min_num;
    den_range = max_den / min_den;
    
    % Apply scaling if dynamic range is too large
    if num_range > 1e8 || den_range > 1e8
        % Calculate overall scaling to reduce dynamic range
        max_coeff = max(max_num, max_den);
        
        if max_coeff > 1e6
            % Scale down if coefficients are too large
            scale = 1e6 / max_coeff;
            num = num * scale;
        elseif max_coeff < 1e-6
            % Scale up if coefficients are too small
            scale = 1e-6 / max_coeff;
            num = num * scale;
        end
    end
end

function coeff = roundToZero(coeff, threshold)
    % Round coefficients below threshold to zero
    
    coeff(abs(coeff) < threshold) = 0;
end

function [num, den] = normalizeCoefficients(num, den)
    % Normalize coefficients to have leading denominator coefficient = 1
    
    % Check for empty or all-zero arrays
    if isempty(den) || all(den == 0)
        return;
    end
    
    % Find leading non-zero coefficient in denominator
    idx = find(abs(den) > eps, 1, 'first');
    
    if ~isempty(idx)
        scale = den(idx);
        den = den / scale;
        num = num / scale;
    end
end

function sys_balanced = getBalancedRealization(sys)
    % GETBALANCEDREALIZATION Get numerically well-conditioned balanced realization
    % of state-space system with enhanced robustness for difficult systems
    %
    % Input:
    %   sys - Original system (state-space or transfer function)
    %
    % Output:
    %   sys_balanced - Balanced realization with improved numerical properties
    
    % Convert to state-space if needed
    if ~isa(sys, 'ss')
        try
            sys_ss = ss(sys);
        catch ME
            warning('State-space conversion failed: %s\nTrying robust conversion.', ME.message);
            sys_ss = getRobustStateSpace(sys);
        end
    else
        sys_ss = sys;
    end
    
    % Try standard balancing
    try
        sys_balanced = balreal(sys_ss);
        % Check if balancing improved conditioning
        if isBetterConditioned(sys_balanced, sys_ss)
            return;
        end
    catch ME
        warning('Standard balancing failed: %s\nTrying alternative methods.', ME.message);
    end
    
    % Try alternative balancing methods
    try
        % Method 1: Schur balancing
        [A, B, C, D] = ssdata(sys_ss);
        [U, T] = schur(A, 'real');
        Tinv = U';
        
        % Transform system
        Ab = Tinv * A * U;
        Bb = Tinv * B;
        Cb = C * U;
        Db = D;
        
        sys_balanced = ss(Ab, Bb, Cb, Db);
        
        % Check if this improved conditioning
        if isBetterConditioned(sys_balanced, sys_ss)
            return;
        end
    catch
        % Continue to next method
    end
    
    % Method 2: Modal form
    try
        sys_modal = ss(sys_ss, 'modal');
        
        % Check if this improved conditioning
        if isBetterConditioned(sys_modal, sys_ss)
            sys_balanced = sys_modal;
            return;
        end
    catch
        % Continue to next method
    end
    
    % Method 3: Apply direct scaling to improve conditioning
    try
        [A, B, C, D] = ssdata(sys_ss);
        n = size(A, 1);
        
        % Create diagonal scaling matrix
        S = eye(n);
        
        for i = 1:n
            % Calculate appropriate scaling for each state
            row_norm = norm(A(i,:));
            col_norm = norm(A(:,i));
            
            if row_norm > 0 && col_norm > 0
                scale = sqrt(col_norm / row_norm);
                
                % Limit scaling to avoid numerical issues
                scale = min(max(scale, 1e-6), 1e6);
                
                S(i,i) = scale;
            end
        end
        
        % Apply scaling transformation
        Sinv = inv(S);
        A_scaled = Sinv * A * S;
        B_scaled = Sinv * B;
        C_scaled = C * S;
        
        sys_balanced = ss(A_scaled, B_scaled, C_scaled, D);
        
        % Check if this improved conditioning
        if isBetterConditioned(sys_balanced, sys_ss)
            return;
        end
    catch
        % Fall back to original system
        sys_balanced = sys_ss;
    end
end

function better = isBetterConditioned(sys1, sys2)
    % Check if sys1 has better numerical conditioning than sys2
    
    [A1, ~, ~, ~] = ssdata(sys1);
    [A2, ~, ~, ~] = ssdata(sys2);
    
    % Calculate condition numbers
    cond1 = cond(A1);
    cond2 = cond(A2);
    
    % Check if condition number improved
    better = cond1 < cond2;
end

function sys_ss = getRobustStateSpace(G)
    % Create a robust state-space representation from transfer function
    % This function handles difficult cases when standard ss() fails
    
    try
        % Get transfer function data
        [num, den] = tfdata(G, 'v');
        
        % First clean and condition the transfer function
        [num_c, den_c] = conditionTransferFunction(num, den);
        
        % Try with improved coefficients
        G_clean = tf(num_c, den_c);
        try
            sys_ss = ss(G_clean);
            return;
        catch
            % Continue with more robust methods
        end
        
        % Try using controllable canonical form directly
        try
            [A, B, C, D] = tf2ss(num_c, den_c);
            sys_ss = ss(A, B, C, D);
            return;
        catch
            % Continue with more robust methods
        end
        
        % Try to decompose into partial fractions and build state-space from parts
        try
            % First analyze poles
            poles = roots(den_c);
            
            % Check for repeated poles
            [unique_poles, ~, multiplicity] = findMultiplePoles(poles);
            
            % Build state-space system using modal form
            n = length(poles);
            A = zeros(n, n);
            B = zeros(n, 1);
            C = zeros(1, n);
            
            idx = 1;
            for i = 1:length(unique_poles)
                p = unique_poles(i);
                m = multiplicity(i);
                
                if m == 1
                    % Simple pole
                    A(idx, idx) = p;
                    B(idx) = 1;
                    
                    % Find residue
                    r = calculateResidue(G_clean, p);
                    C(idx) = r;
                    
                    idx = idx + 1;
                else
                    % Repeated pole - use Jordan block
                    block_idx = idx:(idx+m-1);
                    
                    % Create Jordan block
                    A(block_idx, block_idx) = diag(p*ones(1,m)) + diag(ones(1,m-1), 1);
                    B(idx) = 1;
                    
                    % Calculate residue and derivatives
                    for j = 1:m
                        r = calculateResidue(G_clean, p, j-1);
                        C(idx+j-1) = r;
                    end
                    
                    idx = idx + m;
                end
            end
            
            % Add direct feedthrough if needed
            if length(num_c) >= length(den_c)
                D = num_c(1) / den_c(1);
            else
                D = 0;
            end
            
            sys_ss = ss(A, B, C, D);
            return;
        catch
            % Continue with simplest approach
        end
        
        % Last resort: Simplify the model and try again
        try
            % Reduce model order if needed
            max_order = 10;  % Maximum manageable order
            if length(den_c) > max_order + 1
                warning('Reducing order for numerical stability from %d to %d', length(den_c)-1, max_order);
                
                % Create simplified transfer function
                order = min(length(den_c)-1, max_order);
                [b, a] = prony(impulse(G_clean), length(impulse(G_clean)), order, order);
                
                G_simple = tf(b, a);
                sys_ss = ss(G_simple);
                return;
            end
        catch
            % Ultimate fallback - create trivial system
            warning('All attempts to create state-space model failed. Creating simple approximation.');
            
            % Create simple first-order approximation
            try
                dc = dcgain(G);
                if isnan(dc) || isinf(dc)
                    dc = 1;
                end
            catch
                dc = 1;
            end
            
            sys_ss = ss(-1, 1, dc, 0);
        end
    catch ME
        error('Failed to create robust state-space model: %s', ME.message);
    end
end

function [unique_poles, indices, multiplicity] = findMultiplePoles(poles)
    % Find unique poles and their multiplicity
    
    % Sort poles for consistency
    poles = poles(:);
    
    % Use a tolerance for numerical comparison
    tol = 1e-6;
    
    n = length(poles);
    visited = false(n, 1);
    unique_poles = [];
    indices = {};
    multiplicity = [];
    
    for i = 1:n
        if ~visited(i)
            % Find all occurrences of this pole
            idx = find(abs(poles - poles(i)) < tol);
            visited(idx) = true;
            
            unique_poles = [unique_poles; poles(i)];
            indices{end+1} = idx;
            multiplicity = [multiplicity; length(idx)];
        end
    end
end

function residue = calculateResidue(G, pole, derivative)
    % Calculate the residue of a pole in a transfer function
    % For repeated poles, use the derivative parameter
    
    if nargin < 3
        derivative = 0;
    end
    
    [num, den] = tfdata(G, 'v');
    
    if derivative == 0
        % Standard residue for simple pole
        s = pole;
        
        % Evaluate G(s) * (s - pole)
        n_val = polyval(num, s);
        d_val = polyval(den, s);
        d_der = polyval(polyder(den), s);
        
        residue = n_val / (d_val / (s - pole) + d_der);
    else
        % For higher derivatives of repeated poles
        error('Higher-order residue calculation not implemented yet');
    end
end

function [K_robust, details] = makeControllerRobust(K, details)
    % MAKECONTROLLERROBUST Improve numerical robustness of a controller
    % Applies several techniques to make the controller more robust to
    % numerical issues during implementation
    %
    % Inputs:
    %   K       - Original controller transfer function
    %   details - Design details string
    %
    % Outputs:
    %   K_robust - Numerically robust controller
    %   details  - Updated details with robustness information
    
    details = [details, '\nNUMERICAL ROBUSTNESS IMPROVEMENTS\n'];
    details = [details, '--------------------------------\n'];
    
    % Get controller data
    [num, den] = tfdata(K, 'v');
    
    % 1. Check for numerical conditioning issues
    max_coeff = max(max(abs(num)), max(abs(den)));
    min_num = min(abs(num(abs(num) > eps)));
    min_den = min(abs(den(abs(den) > eps)));
    
    if isempty(min_num), min_num = 1; end
    if isempty(min_den), min_den = 1; end
    
    min_coeff = min(min_num, min_den);
    
    condition_number = max_coeff / min_coeff;
    
    details = [details, sprintf('Original controller condition number: %.2e\n', condition_number)];
    
    if condition_number > 1e8
        details = [details, 'Controller has poor numerical conditioning. Applying improvements.\n'];
        
        % 2. Apply numerical conditioning
        [num_c, den_c] = conditionTransferFunction(num, den);
        
        % 3. Check for pole-zero cancellations
        [num_c, den_c, cancellations] = removePoleZeroCancellations(num_c, den_c);
        
        if cancellations > 0
            details = [details, sprintf('Removed %d approximate pole-zero cancellations.\n', cancellations)];
        end
        
        % 4. Create the robust controller
        K_robust = tf(num_c, den_c);
        
        % 5. Verify improvement
        max_coeff_new = max(max(abs(num_c)), max(abs(den_c)));
        min_num_new = min(abs(num_c(abs(num_c) > eps)));
        min_den_new = min(abs(den_c(abs(den_c) > eps)));
        
        if isempty(min_num_new), min_num_new = 1; end
        if isempty(min_den_new), min_den_new = 1; end
        
        min_coeff_new = min(min_num_new, min_den_new);
        
        condition_number_new = max_coeff_new / min_coeff_new;
        
        details = [details, sprintf('Improved controller condition number: %.2e\n', condition_number_new)];
        details = [details, sprintf('Improvement factor: %.2f\n', condition_number / condition_number_new)];
    else
        details = [details, 'Controller has good numerical conditioning. No changes needed.\n'];
        K_robust = K;
    }
    
    % 6. Check for controller stability
    p = pole(K_robust);
    unstable_poles = sum(real(p) >= 0);
    
    if unstable_poles > 0
        details = [details, sprintf('WARNING: Controller has %d pole(s) in RHP.\n', unstable_poles)];
        details = [details, 'This is normal for some design methods but may cause implementation issues.\n'];
        details = [details, 'Consider using state-variable implementation for this controller.\n'];
    else
        details = [details, 'Controller is stable (all poles in LHP). Good for implementation.\n'];
    }
    
    % 7. Check for high-frequency characteristics
    try
        [mag, ~, w] = bode(K_robust, {1e0, 1e4});
        mag = squeeze(mag);
        hf_gain = mag(end);
        
        if hf_gain > 100
            details = [details, sprintf('WARNING: Controller has high gain (%.1f dB) at high frequencies.\n', 20*log10(hf_gain))];
            details = [details, 'This may amplify measurement noise. Consider additional filtering if needed.\n'];
        end
    catch
        % Skip high-frequency check if it fails
    end
    
    % 8. Recommend discretization method if appropriate
    details = [details, '\nIMPLEMENTATION RECOMMENDATIONS:\n'];
    
    if unstable_poles > 0
        details = [details, '- Use state-space implementation for this controller\n'];
        details = [details, '- Use bilinear (Tustin) method with prewarping if discretizing\n'];
        details = [details, '- Verify stability after discretization\n'];
    else
        details = [details, '- Standard implementation methods are suitable\n'];
        details = [details, '- For discretization, recommend bilinear (Tustin) method\n'];
    end
    
    if any(abs(imag(p)) > 10 * abs(real(p)))
        details = [details, '- Controller has lightly damped poles. Use high sampling rate if discretizing\n'];
    end
    
    % 9. Add details about minreal tolerance used
    details = [details, sprintf('- Pole-zero cancellation tolerance: %.2e\n', 1e-4)];
end

function [num, den, cancellations] = removePoleZeroCancellations(num, den)
    % Remove approximate pole-zero cancellations for better numerical conditioning
    
    % Initial count
    cancellations = 0;
    
    % Calculate poles and zeros
    z = roots(num);
    p = roots(den);
    
    % Skip if empty
    if isempty(z) || isempty(p)
        return;
    end
    
    % Use tolerance for matching
    tol = 1e-4;
    
    % Track which poles and zeros have been cancelled
    z_cancelled = false(size(z));
    p_cancelled = false(size(p));
    
    % Find pole-zero pairs that are close
    for i = 1:length(z)
        if z_cancelled(i)
            continue;
        end
        
        for j = 1:length(p)
            if p_cancelled(j)
                continue;
            end
            
            % Check if this is a cancellation
            if abs(z(i) - p(j)) < tol
                z_cancelled(i) = true;
                p_cancelled(j) = true;
                cancellations = cancellations + 1;
                break;
            end
        end
    end
    
    % Create new polynomials without the cancelled terms
    z_remain = z(~z_cancelled);
    p_remain = p(~p_cancelled);
    
    % Create new numerator and denominator
    if isempty(z_remain)
        num_new = 1;  % Constant term
    else
        num_new = poly(z_remain);
    end
    
    if isempty(p_remain)
        den_new = 1;  % Should not happen normally
    else
        den_new = poly(p_remain);
    end
    
    % Preserve DC gain
    dc_orig = polyval(num, 0) / polyval(den, 0);
    dc_new = polyval(num_new, 0) / polyval(den_new, 0);
    
    % Adjust gain
    num = num_new * (dc_orig / dc_new);
    den = den_new;
end