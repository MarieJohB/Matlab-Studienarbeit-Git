function [K, details] = designMuSynthesis(G, structure, options, plantInfo)
% DESIGNMUSYNTHESIS Controller design using robust µ-synthesis approach
%
% Implements an advanced µ-synthesis controller design for systems with uncertainty.
% This method has been enhanced to better handle high-order unstable systems and
% to integrate state-space model information when available.
%
% Inputs:
%   G        - Plant transfer function or state-space model
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters
%   plantInfo - Structure with plant analysis information
%
% Outputs:
%   K       - Designed controller as transfer function
%   details - Text description with design details

    % Start with detailed information about the method
    details = 'Robust µ-Synthesis Design Method\n';
    details = [details, '------------------------------\n'];
    details = [details, sprintf('Plant Analysis: %s\n', getPlantInfoString(plantInfo))];

    % Default values
    if ~isfield(options, 'bandwidth')
        options.bandwidth = 1;
    end
    
    if ~isfield(options, 'damping')
        options.damping = 0.8;
    end
    
    if ~isfield(options, 'epsilon')
        options.epsilon = 0.1;
    end
    
    if ~isfield(options, 'robustness')
        options.robustness = 'Medium';
    end
    
    if ~isfield(options, 'uncertainty')
        options.uncertainty = 20;  % Default 20% uncertainty
    end
    
    % Extract key parameters
    omega = options.bandwidth;
    zeta = options.damping;
    epsilon = options.epsilon;
    robustness = options.robustness;
    uncertaintyPercent = options.uncertainty;
    
    details = [details, sprintf('Desired bandwidth: %.4f rad/s\n', omega)];
    details = [details, sprintf('Desired damping ratio: %.4f\n', zeta)];
    details = [details, sprintf('Derivative filter coefficient: %.4f\n', epsilon)];
    details = [details, sprintf('Robustness level: %s\n', robustness)];
    details = [details, sprintf('Uncertainty percentage: %.1f%%\n', uncertaintyPercent)];
    
    % Check if input is state-space or transfer function
    isStateSpace = isa(G, 'ss');
    
    % If state-space model is in the options, use it directly
    if ~isStateSpace && isfield(options, 'stateSpace')
        G_ss = options.stateSpace;
        isStateSpace = true;
        details = [details, 'Using state-space model from options.\n'];
    elseif isStateSpace
        G_ss = G;
        details = [details, 'Using provided state-space model directly.\n'];
    else
        % Try to convert to state-space for more robust handling
        try
            G_ss = ss(G);
            isStateSpace = true;
            details = [details, 'Successfully converted to state-space representation.\n'];
        catch ME
            details = [details, sprintf('Could not convert to state-space: %s\n', ME.message)];
            details = [details, 'Using transfer function approach.\n'];
            G_ss = [];
        end
    end
    
    % Check if Robust Control Toolbox is available
    hasRobustToolbox = exist('musyn', 'file') == 2;
    
    if ~hasRobustToolbox
        details = [details, 'Robust Control Toolbox not available. Using D-K iteration approximation.\n'];
    else
        details = [details, 'Robust Control Toolbox available. Using direct µ-synthesis.\n'];
    end
    
    % Step 1: Create uncertain plant model
    details = [details, '\nSTEP 1: Creating uncertain plant model\n'];
    details = [details, '-----------------------------------\n'];
    
    % Special handling for highly unstable systems
    isHighlyUnstable = false;
    if plantInfo.isUnstable
        p = plantInfo.poles;
        unstable_poles = p(real(p) > 0);
        
        if max(real(unstable_poles)) > 5 || length(unstable_poles) > 1
            isHighlyUnstable = true;
            details = [details, 'System is highly unstable. Using specialized approach.\n'];
        end
    end
    
    try
        % Define uncertainty based on plant characteristics
        if hasRobustToolbox
            details = [details, 'Creating structured uncertainty model...\n'];
            
            % Create nominal plant model
            if isStateSpace
                G_nom = G_ss;
            else
                G_nom = G;
            end
            
            % Define uncertainty level based on robustness setting
            switch robustness
                case 'Low'
                    delta_gain = uncertaintyPercent / 100;  % Less conservative
                case 'High'
                    delta_gain = 2 * uncertaintyPercent / 100;  % More conservative
                otherwise  % Medium
                    delta_gain = 1.5 * uncertaintyPercent / 100;
            end
            
            % For highly unstable systems, increase uncertainty level
            if isHighlyUnstable
                delta_gain = delta_gain * 1.5;
                details = [details, 'Increased uncertainty level for highly unstable system.\n'];
            end
            
            % Create input and output uncertainty models
            if plantInfo.isUnstable
                details = [details, 'Using multiplicative input and output uncertainty models for unstable plant.\n'];
                
                % For unstable plants, use more sophisticated uncertainty modeling
                Wout = makeweight(delta_gain/3, omega*1.5, delta_gain);
                Win = makeweight(delta_gain/3, omega*1.5, delta_gain);
                
                % Add more uncertainty for non-minimum phase plants
                if plantInfo.hasRHPZeros
                    Wout = Wout * 1.2;
                    details = [details, 'Increased uncertainty for non-minimum phase behavior.\n'];
                end
                
                % Create uncertain plant
                Delta_in = ultidyn('Delta_in', [1 1], 'Bound', 1);
                Delta_out = ultidyn('Delta_out', [1 1], 'Bound', 1);
                
                G_unc = G_nom * (1 + Win*Delta_in) * (1 + Wout*Delta_out);
                
            else
                % For stable plants, simpler uncertainty model is sufficient
                details = [details, 'Using multiplicative input uncertainty model for stable plant.\n'];
                
                % Create input uncertainty weight
                Win = makeweight(delta_gain/4, omega, delta_gain/2);
                
                % Add more uncertainty for non-minimum phase plants
                if plantInfo.hasRHPZeros
                    Win = Win * 1.2;
                    details = [details, 'Increased uncertainty for non-minimum phase behavior.\n'];
                end
                
                % Create uncertain plant
                Delta_in = ultidyn('Delta_in', [1 1], 'Bound', 1);
                G_unc = G_nom * (1 + Win*Delta_in);
            end
        else
            % Without Robust Control Toolbox, we'll use manual uncertainty representation
            details = [details, 'Creating manual uncertainty representation for D-K iteration...\n'];
            
            % Define uncertainty weights based on plant characteristics
            delta_gain = uncertaintyPercent / 100;
            
            % Adjust based on robustness setting
            switch robustness
                case 'Low'
                    delta_gain = delta_gain * 0.8;  % Less conservative
                case 'High'
                    delta_gain = delta_gain * 1.5;  % More conservative
            end
            
            % For highly unstable systems, increase uncertainty level
            if isHighlyUnstable
                delta_gain = delta_gain * 1.5;
                details = [details, 'Increased uncertainty level for highly unstable system.\n'];
            end
            
            % Create appropriate frequency-dependent weights
            if plantInfo.isUnstable
                % For unstable plants, more sophisticated weighting
                Wout = tf([delta_gain/3, delta_gain*omega*1.5], [1, omega*1.5/delta_gain/3]);
                Win = tf([delta_gain/3, delta_gain*omega*1.5], [1, omega*1.5/delta_gain/3]);
                
                % Add more uncertainty for non-minimum phase plants
                if plantInfo.hasRHPZeros
                    Wout = Wout * 1.2;
                    details = [details, 'Increased uncertainty for non-minimum phase behavior.\n'];
                end
            else
                % For stable plants
                Win = tf([delta_gain/4, delta_gain*omega/2], [1, omega/delta_gain/4]);
                Wout = tf(0.01, 1);  % Minimal output uncertainty for stable plants
                
                % Add more uncertainty for non-minimum phase plants
                if plantInfo.hasRHPZeros
                    Win = Win * 1.2;
                    details = [details, 'Increased uncertainty for non-minimum phase behavior.\n'];
                end
            end
            
            if isStateSpace
                G_unc = G_ss;  % Use state-space model
            else
                G_unc = G;  % Use transfer function model
            end
        end
        
    catch ME
        details = [details, sprintf('Error creating uncertainty model: %s\n', ME.message)];
        details = [details, 'Using simplified uncertainty model...\n'];
        
        % Create basic uncertainty weights
        delta_gain = uncertaintyPercent / 100;
        
        % Simple first-order input uncertainty weight
        Win = tf([delta_gain, delta_gain*omega], [1, omega/delta_gain]);
        
        % Simple scalar output uncertainty
        Wout = tf(delta_gain/10, 1);
        
        if isStateSpace
            G_unc = G_ss;
        else
            G_unc = G;
        end
    end
    
    % Step 2: Create generalized plant P for synthesis
    details = [details, '\nSTEP 2: Creating generalized plant P for synthesis\n'];
    details = [details, '---------------------------------------------\n'];
    
    try
        % Create performance weights based on desired specification
        
        % Weight on sensitivity function (S = 1/(1+GK))
        % Lower S means better tracking and disturbance rejection
        switch robustness
            case 'Low'
                % Lower robustness: emphasize performance
                M_s = 2.0;  % Higher peak sensitivity allowed
                omega_b = omega / 10;  % Bandwidth for S
                A_s = 0.01; % Steady-state error
                
                % Create weight with lower DC gain but wider bandwidth
                Ws = tf([1/sqrt(M_s), omega_b], [1, omega_b*A_s]);
                
            case 'High'
                % Higher robustness: more conservative
                M_s = 1.2;  % Lower peak sensitivity
                omega_b = omega / 20;  % Reduced bandwidth
                A_s = 0.1;  % Larger allowed steady-state error
                
                % Create weight with higher DC gain but narrower bandwidth
                Ws = tf([1/sqrt(M_s), omega_b], [1, omega_b*A_s]);
                
                % Add extra roll-off for highly robust controllers
                Ws = Ws * tf([1, omega*2], [1, omega*10]);
                
            otherwise  % Medium
                % Balanced approach
                M_s = 1.5;  % Moderate peak sensitivity
                omega_b = omega / 15;  % Moderate bandwidth
                A_s = 0.05;  % Moderate steady-state error
                
                Ws = tf([1/sqrt(M_s), omega_b], [1, omega_b*A_s]);
        end
        
        % Adjust Ws for plants with integrator
        if plantInfo.hasIntegrator
            % Reduce the DC gain requirement since plant already has integrator
            Ws = Ws * tf([1, omega_b/10], [1, omega_b/100]);
            details = [details, 'Adjusted sensitivity weight for plant with integrator.\n'];
        end
        
        % Weight on control sensitivity function (KS)
        % Limits control effort and improves robustness to high-frequency unmodeled dynamics
        Wu = tf(1, 100);  % Base low-effort weight
        
        % Adjust based on plant characteristics
        if plantInfo.isHighOrder
            % For high-order plants, limit high-frequency control action more
            Wu = Wu * tf([1, omega*2], [1, omega*20]);
        end
        
        if plantInfo.isUnstable
            % For unstable plants, allow more control effort at critical frequencies
            p = plantInfo.poles;
            unstable_poles = p(real(p) > 0);
            
            if ~isempty(unstable_poles)
                fastest_unstable = max(abs(unstable_poles));
                Wu = Wu * tf([1, fastest_unstable*2], [1, fastest_unstable/2]);
                details = [details, sprintf('Adjusted control effort weight around unstable pole frequencies (%.3f rad/s).\n', fastest_unstable)];
            end
            
            % For highly unstable systems, allow even more control effort
            if isHighlyUnstable
                Wu = Wu / 2;
                details = [details, 'Allowed more control effort for highly unstable system.\n'];
            end
        end
        
        % Weight on complementary sensitivity function (T = GK/(1+GK))
        % Ensures robustness to unmodeled dynamics and noise rejection
        switch robustness
            case 'Low'
                % Less high-frequency roll-off
                Wt = tf([1, omega], [0.1, omega*5]);
                
            case 'High'
                % More high-frequency roll-off and lower bandwidth
                Wt = tf([1, omega/2], [0.01, omega*10]);
                
                % Additional roll-off for more robust designs
                Wt = Wt * tf([1, omega*2], [1, omega*10]);
                
            otherwise  % Medium
                Wt = tf([1, omega/1.5], [0.05, omega*8]);
        end
        
        % Adjust Wt for non-minimum phase plants
        if plantInfo.hasRHPZeros
            z = plantInfo.zeros;
            rhp_zeros = z(real(z) > 0);
            
            if ~isempty(rhp_zeros)
                min_rhp_zero = min(real(rhp_zeros));
                
                % Reduce bandwidth required close to RHP zeros
                Wt = Wt * tf([1, min_rhp_zero/2], [1, min_rhp_zero*2]);
                details = [details, sprintf('Adjusted complementary sensitivity weight due to RHP zero at %.3f.\n', min_rhp_zero)];
            end
        end
        
        % Create augmented plant for synthesis
        details = [details, '\nCreating augmented plant with performance weights:\n'];
        details = [details, sprintf('Tracking/Sensitivity weight Ws: %s\n', char(Ws))];
        details = [details, sprintf('Control effort weight Wu: %s\n', char(Wu))];
        details = [details, sprintf('Robustness weight Wt: %s\n', char(Wt))];
        
        if hasRobustToolbox
            % Create generalized plant for mu-synthesis with special handling for unstable plants
            if plantInfo.isUnstable && plantInfo.hasRHPZeros
                details = [details, 'Using specialized interconnection for unstable non-minimum phase plant.\n'];
                
                systemnames = 'G_unc Ws Wu Wt';
                inputvar = '[r; d; n; u]';
                outputvar = '[Ws; Wu; Wt; r-G_unc]';
                input_to_G_unc = '[u]';
                input_to_Ws = '[r-G_unc]';
                input_to_Wu = '[u]';
                input_to_Wt = '[G_unc]';
                
                try
                    P = sysic;
                catch ME
                    details = [details, sprintf('sysic failed: %s\n', ME.message)];
                    
                    if isStateSpace
                        P = augw(G_ss, Ws, Wu, Wt);
                    else
                        P = augw(G, Ws, Wu, Wt);
                    end
                    
                    details = [details, 'Using augw as fallback.\n'];
                end
            else
                % Standard mixed-sensitivity setup
                if isStateSpace
                    P = augw(G_ss, Ws, Wu, Wt);
                else
                    P = augw(G_unc, Ws, Wu, Wt);
                end
            end
        else
            % Manual creation of generalized plant without Robust Control Toolbox
            details = [details, 'Creating manual generalized plant for D-K iteration.\n'];
            
            if isStateSpace
                P = augw(G_ss, Ws, Wu, Wt);
            else
                P = augw(G, Ws, Wu, Wt);
            end
        end
    catch ME
        details = [details, sprintf('Error creating generalized plant: %s\n', ME.message)];
        details = [details, 'Using simplified generalized plant...\n'];
        
        % Create a simplified generalized plant for synthesis
        % Simple mixed-sensitivity setup
        Ws = tf([1/1.5, omega/10], [1, omega/1000]);  % Tracking
        Wu = tf(0.1, 1);  % Control effort
        Wt = tf([1, omega], [0.1, omega*5]);  % Complementary sensitivity
        
        % Use loopsyn for simple plant creation
        if isStateSpace
            P = augw(G_ss, Ws, Wu, Wt);
        else
            P = augw(G, Ws, Wu, Wt);
        end
    end
    
    % Step 3: Perform mu-synthesis or D-K iteration
    details = [details, '\nSTEP 3: Performing controller synthesis\n'];
    details = [details, '-------------------------------------\n'];
    
    try
        if hasRobustToolbox
            % Use direct mu-synthesis
            details = [details, 'Using musyn for direct µ-synthesis...\n'];
            
            % Options for musyn
            opt = musynOptions;
            opt.MaxIter = 5;
            
            % For highly unstable systems, adjust options
            if isHighlyUnstable
                opt.MaxIter = 8;  % More iterations
                details = [details, 'Increased maximum iterations for highly unstable system.\n'];
            end
            
            % Perform mu-synthesis
            [K_mu, CL, mubnds] = musyn(P, 4, 1, opt);
            
            details = [details, sprintf('µ-synthesis completed with %d D-K iterations.\n', length(mubnds))];
            details = [details, sprintf('Final µ-bound: %.4f\n', mubnds(end))];
            
            % Use the resulting controller
            K_robust = K_mu;
        else
            % Implement manual D-K iteration approximation
            details = [details, 'Performing manual D-K iteration approximation...\n'];
            
            % Step 1: Initial H-infinity design
            try
                [K1, CL, gamma] = hinfsyn(P, 1, 1);
                details = [details, sprintf('Initial H-infinity design complete with gamma = %.4f\n', gamma)];
            catch ME
                details = [details, sprintf('Initial H-infinity design failed: %s\n', ME.message)];
                details = [details, 'Using loopsyn as fallback for initial design.\n'];
                
                if isStateSpace
                    K1 = loopsyn(G_ss, Ws);
                else
                    K1 = loopsyn(G, Ws);
                end
                
                gamma = NaN;
            end
            
            % Step 2: Analyze and adjust weights
            % We'll simulate D-K iteration by adjusting weights and redesigning
            try
                % Extract frequency response of initial controller
                w = logspace(-2, log10(omega*20), 100);
                
                try
                    if exist('CL', 'var')
                        mag = sigma(CL, w);
                    else
                        if isStateSpace
                            mag = sigma(feedback(G_ss*K1, 1), w);
                        else
                            mag = sigma(feedback(G*K1, 1), w);
                        end
                    end
                    
                    mag = squeeze(mag);
                    
                    % Find peak frequency
                    [peak_mag, idx] = max(mag);
                    peak_freq = w(idx);
                    
                    details = [details, sprintf('Peak magnitude %.4f at frequency %.4f rad/s\n', peak_mag, peak_freq)];
                    
                    % Adjust weights to focus on critical frequencies
                    Ws_adj = Ws * tf([1, peak_freq*0.5], [1, peak_freq*2]);
                    Wt_adj = Wt * tf([1, peak_freq*0.5], [1, peak_freq*2]);
                catch
                    % If analysis fails, use simpler adjustment
                    Ws_adj = Ws * 1.2;
                    Wt_adj = Wt * 1.2;
                    details = [details, 'Using simple weight scaling for second iteration.\n'];
                end
                
                % Create adjusted plant
                if isStateSpace
                    P_adj = augw(G_ss, Ws_adj, Wu, Wt_adj);
                else
                    P_adj = augw(G, Ws_adj, Wu, Wt_adj);
                end
                
                % Second iteration
                try
                    [K2, CL2, gamma2] = hinfsyn(P_adj, 1, 1);
                    details = [details, sprintf('Second iteration complete with gamma = %.4f\n', gamma2)];
                catch
                    details = [details, 'Second iteration failed. Using loopsyn as fallback.\n'];
                    
                    if isStateSpace
                        K2 = loopsyn(G_ss, Ws_adj);
                    else
                        K2 = loopsyn(G, Ws_adj);
                    end
                    
                    gamma2 = NaN;
                end
                
                % Third iteration with further adjustment
                if isnan(gamma) || isnan(gamma2) || gamma2 < gamma
                    % If improvement or one failed, continue in same direction
                    Ws_adj2 = Ws_adj * tf([1, omega*0.2], [1, omega*5]);
                    Wt_adj2 = Wt_adj * tf([1, omega*0.2], [1, omega*5]);
                else
                    % If no improvement, try different approach
                    Ws_adj2 = Ws * tf([1, omega*0.5], [1, omega*5]);
                    Wt_adj2 = Wt * tf([1, omega*2], [1, omega*0.5]);
                end
                
                if isStateSpace
                    P_adj2 = augw(G_ss, Ws_adj2, Wu, Wt_adj2);
                else
                    P_adj2 = augw(G, Ws_adj2, Wu, Wt_adj2);
                end
                
                try
                    [K3, CL3, gamma3] = hinfsyn(P_adj2, 1, 1);
                    details = [details, sprintf('Third iteration complete with gamma = %.4f\n', gamma3)];
                catch
                    details = [details, 'Third iteration failed. Using loopsyn as fallback.\n'];
                    
                    if isStateSpace
                        K3 = loopsyn(G_ss, Ws_adj2);
                    else
                        K3 = loopsyn(G, Ws_adj2);
                    end
                    
                    gamma3 = NaN;
                end
                
                % Select best controller based on performance
                % Check closed-loop stability first
                controllers = {K1, K2, K3};
                stabilities = zeros(1, 3);
                scores = zeros(1, 3);
                
                for i = 1:3
                    try
                        if isStateSpace
                            cl = feedback(G_ss*controllers{i}, 1);
                        else
                            cl = feedback(G*controllers{i}, 1);
                        end
                        
                        % Check if stable
                        if all(real(pole(cl)) < 0)
                            stabilities(i) = 1;
                            
                            % Calculate a score based on frequency response
                            try
                                [mag_s, ~] = bode(cl, w);
                                mag_s = squeeze(mag_s);
                                
                                % Score based on peak and low-frequency performance
                                peak = max(mag_s);
                                low_freq = mag_s(1);
                                
                                scores(i) = 1/(peak * low_freq);
                            catch
                                % If frequency analysis fails, give a positive but lower score
                                scores(i) = 1;
                            end
                        end
                    catch
                        % If analysis fails, mark as unstable
                        stabilities(i) = 0;
                        scores(i) = 0;
                    end
                end
                
                % Choose the best stable controller
                if any(stabilities)
                    % Get the best among stable controllers
                    [~, best_idx] = max(scores .* stabilities);
                    K_robust = controllers{best_idx};
                    details = [details, sprintf('Using controller from iteration %d based on stability and performance.\n', best_idx)];
                else
                    % If none are stable, return to pre-stabilization approach
                    details = [details, 'No stable controllers found. Using pre-stabilization approach.\n'];
                    [K_robust, ~] = preStabilize(G, plantInfo);
                end
            catch ME2
                details = [details, sprintf('D-K iteration process failed: %s\n', ME2.message)];
                details = [details, 'Using initial H-infinity controller as fallback.\n'];
                
                K_robust = K1;
            end
        end
        
        % Get order of the synthesized controller
        [A_K, B_K, C_K, D_K] = ssdata(K_robust);
        controller_order = size(A_K, 1);
        details = [details, sprintf('Controller order: %d\n', controller_order)];
        
        % Step 4: Controller order reduction and conversion to desired structure
        details = [details, '\nSTEP 4: Controller order reduction and structure conversion\n'];
        details = [details, '----------------------------------------------------\n'];
        
        % First, try to reduce controller order
        try
            % Balance the controller for better numerical properties
            K_bal = balreal(K_robust);
            
            % Get Hankel singular values
            hsv = hsvd(K_bal);
            
            % Determine how many states to keep based on Hankel singular values
            if length(hsv) > 4
                % Find threshold where HSVs drop significantly
                hsv_ratios = hsv(1:end-1) ./ hsv(2:end);
                [~, idx] = max(hsv_ratios);
                
                % Keep at least 2 states, but no more than 4 for practical controllers
                keep_states = min(max(idx, 2), 4);
                
                % Reduce order
                K_red = reduce(K_bal, keep_states);
                
                details = [details, sprintf('Reduced controller order from %d to %d\n', controller_order, keep_states)];
                
                K_robust = K_red;
            else
                details = [details, 'Controller already has low order. No reduction performed.\n'];
            end
            
            % Verify stability after reduction
            k_poles = pole(K_robust);
            if any(real(k_poles) > 0)
                details = [details, 'Warning: Reduced controller has unstable poles. Using original controller.\n'];
                K_robust = K_bal;
            end
        catch ME
            details = [details, sprintf('Order reduction failed: %s\n', ME.message)];
            details = [details, 'Continuing with original controller.\n'];
        end
        
    catch ME
        details = [details, sprintf('Error in controller synthesis: %s\n', ME.message)];
        details = [details, 'Using fallback control design approach...\n'];
        
        % For unstable plants, try pre-stabilization
        if plantInfo.isUnstable
            details = [details, 'Using pre-stabilization approach for unstable plant.\n'];
            
            [K_robust, prestab_details] = preStabilize(G, plantInfo);
            details = [details, '-- Pre-stabilization details --\n'];
            details = [details, prestab_details];
        else
            % For stable plants, use loop-shaping approach
            
            % Adjust parameters based on robustness setting
            switch robustness
                case 'Low'
                    phaseMargin = 35;  % Lower phase margin for less robustness
                case 'High'
                    phaseMargin = 60;  % Higher phase margin for more robustness
                otherwise
                    phaseMargin = 45;  % Default phase margin
            end
            
            % Use loop-shaping as fallback
            [K_robust, loop_details] = designLoopShaping(G, structure, phaseMargin, omega, epsilon, plantInfo);
            details = [details, '-- Loop-shaping details --\n'];
            details = [details, loop_details];
        end
    end
    
    % Convert to the desired controller structure
    details = [details, sprintf('\nConverting µ-synthesis controller to %s structure...\n', structure)];
    
    % Verify controller exists and is valid
    if ~exist('K_robust', 'var') || isempty(K_robust)
        details = [details, 'Error: No valid controller was generated.\n'];
        
        % Create a default controller based on structure
        switch structure
            case 'P'
                K = tf(0.5, 1);
            case 'PI'
                K = tf([1, 0.1], [1, 0]);
            case 'PD'
                K = tf([0.1, 1], [0.01, 1]);
            case 'PID'
                K = tf([0.1, 1, 0.1], [0.01, 1, 0]);
            otherwise
                K = tf(0.5, 1);
        end
        
        details = [details, 'Using default controller due to synthesis failure.\n'];
        return;
    end
    
    % Extract controller parameters based on structure
    try
        % Use frequency response for more robust parameter extraction
        w = logspace(-3, log10(omega*10), 200);
        [mag, phase] = bode(K_robust, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        switch structure
            case 'P'
                % Extract proportional gain (at middle frequencies)
                mid_idx = ceil(length(w)/2);
                Kp = mag(mid_idx);
                
                % Limit gain to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                
                K = tf(Kp, 1);
                details = [details, sprintf('Extracted P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                % Extract PI parameters
                
                % Check for integral action (phase approaching -90° at low frequencies)
                has_integrator = (phase(1) < -45);
                
                if has_integrator
                    details = [details, 'Detected integral action in synthesized controller.\n'];
                    
                    % Find frequency where phase crosses -45 degrees
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    % Extract proportional gain at mid-frequencies
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Calculate integral gain
                    Ki = Kp * w_i / 5;
                else
                    details = [details, 'No clear integral action detected. Adding appropriate integral term.\n'];
                    
                    % Extract proportional gain
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Add moderate integral action
                    Ki = Kp * omega / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                        details = [details, 'Reduced integral gain due to plant integrator.\n'];
                    end
                end
                
                % Special handling for unstable systems
                if plantInfo.isUnstable
                    Ki = Ki / 5;
                    details = [details, 'Reduced integral gain for unstable plant.\n'];
                    
                    if isHighlyUnstable
                        Ki = Ki / 5;
                        details = [details, 'Further reduced integral gain for highly unstable plant.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Ki = min(max(abs(Ki), 0.01), 50);
                
                K = tf([Kp, Ki], [1, 0]);
                details = [details, sprintf('Extracted PI controller: Kp = %.4f, Ki = %.4f\n', Kp, Ki)];
                
            case 'PD'
                % Extract PD parameters
                
                % Check for derivative action (phase > 0)
                has_derivative = any(phase > 10);
                
                if has_derivative
                    details = [details, 'Detected derivative action in synthesized controller.\n'];
                    
                    % Find where phase is maximum (derivative action)
                    [~, idx_max] = max(phase);
                    w_d = w(idx_max);
                    
                    % Extract proportional gain at mid-frequencies
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Calculate derivative gain
                    Kd = Kp / w_d;
                else
                    details = [details, 'No strong derivative action detected. Adding appropriate derivative term.\n'];
                    
                    % Extract proportional gain
                    idx_mid = ceil(length(w)/2);
                    Kp = mag(idx_mid);
                    
                    % Add derivative action based on bandwidth
                    Kd = Kp / omega;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                        details = [details, 'Reduced derivative action due to RHP zeros.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Kd = min(max(abs(Kd), 0.01), 50);
                
                % Add filtering for derivative term
                Td = Kd / Kp;
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('Extracted PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
                
            case 'PID'
                % Extract PID parameters
                
                % Check for both integral and derivative actions
                has_integrator = (phase(1) < -45);
                has_derivative = any(phase > 10);
                
                if has_integrator && has_derivative
                    details = [details, 'Detected both integral and derivative action in synthesized controller.\n'];
                    
                    % Find frequencies for integral and derivative components
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    [~, idx_max] = max(phase);
                    w_d = w(idx_max);
                    
                    % Extract parameters
                    Kp = mag(ceil(length(w)/2));  % Mid-frequency gain
                    Ki = Kp * w_i / 5;
                    Kd = Kp / w_d;
                else
                    details = [details, 'Full PID behavior not detected. Creating balanced PID controller.\n'];
                    
                    % Extract proportional gain
                    Kp = mag(ceil(length(w)/2));
                    
                    % Add integral and derivative actions based on bandwidth
                    Ki = Kp * omega / 10;
                    Kd = Kp / omega;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                        details = [details, 'Reduced integral action due to plant integrator.\n'];
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                        details = [details, 'Reduced derivative action due to RHP zeros.\n'];
                    end
                end
                
                % Special handling for unstable systems
                if plantInfo.isUnstable
                    Ki = Ki / 5;
                    details = [details, 'Reduced integral action for unstable plant.\n'];
                    
                    if isHighlyUnstable
                        Ki = Ki / 5;
                        details = [details, 'Further reduced integral gain for highly unstable plant.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(abs(Kp), 0.1), 100);
                Ki = min(max(abs(Ki), 0.01), 50);
                Kd = min(max(abs(Kd), 0.01), 50);
                
                % Create PID controller with derivative filter
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                details = [details, sprintf('Extracted PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', Kp, Ki, Kd)];
                
            otherwise
                error('Unsupported controller structure');
        end
    catch ME
        details = [details, sprintf('Error extracting controller parameters: %s\n', ME.message)];
        details = [details, 'Using simplified parameter extraction...\n'];
        
        % Fallback extraction using DC gain and order
        try
            K_gain = abs(dcgain(K_robust));
            if isnan(K_gain) || isinf(K_gain)
                % Try gain at bandwidth frequency
                K_gain = abs(evalfr(K_robust, 1j*omega));
                
                if isnan(K_gain) || isinf(K_gain)
                    K_gain = 1.0;  % Default
                end
            end
            
            % Limit gain to reasonable values
            K_gain = min(max(K_gain, 0.1), 10);
            
            % Create controller based on structure
            switch structure
                case 'P'
                    K = tf(K_gain, 1);
                    details = [details, sprintf('Created basic P controller: Kp = %.4f\n', K_gain)];
                    
                case 'PI'
                    Ki = K_gain * omega / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                    end
                    
                    if plantInfo.isUnstable
                        Ki = Ki / 5;
                        
                        if isHighlyUnstable
                            Ki = Ki / 5;
                        end
                    end
                    
                    K = tf([K_gain, Ki], [1, 0]);
                    details = [details, sprintf('Created basic PI controller: Kp = %.4f, Ki = %.4f\n', K_gain, Ki)];
                    
                case 'PD'
                    Kd = K_gain / omega;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                    end
                    
                    Td = Kd / K_gain;
                    K = tf([Kd, K_gain], [epsilon*Td, 1]);
                    details = [details, sprintf('Created basic PD controller: Kp = %.4f, Kd = %.4f\n', K_gain, Kd)];
                    
                case 'PID'
                    Ki = K_gain * omega / 10;
                    Kd = K_gain / omega;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                    end
                    
                    if plantInfo.isUnstable
                        Ki = Ki / 5;
                        
                        if isHighlyUnstable
                            Ki = Ki / 5;
                        end
                    end
                    
                    K = tf([Kd, K_gain, Ki], [epsilon*Kd, 1, 0]);
                    details = [details, sprintf('Created basic PID controller: Kp = %.4f, Ki = %.4f, Kd = %.4f\n', K_gain, Ki, Kd)];
                    
                otherwise
                    K = tf(K_gain, 1);
                    details = [details, 'Using fallback P controller.\n'];
            end
        catch
            % If everything fails, create a default controller
            details = [details, 'Creating default safe controller.\n'];
            
            switch structure
                case 'P'
                    K = tf(1, 1);
                case 'PI'
                    if plantInfo.isUnstable
                        K = tf([1, 0.01], [1, 0]);
                    else
                        K = tf([1, 0.1], [1, 0]);
                    end
                case 'PD'
                    K = tf([0.1, 1], [epsilon*0.1, 1]);
                case 'PID'
                    if plantInfo.isUnstable
                        K = tf([0.1, 1, 0.01], [epsilon*0.1, 1, 0]);
                    else
                        K = tf([0.1, 1, 0.1], [epsilon*0.1, 1, 0]);
                    end
                otherwise
                    K = tf(1, 1);
            end
            
            details = [details, 'Using default controller parameters as extraction failed.\n'];
        end
    end
    
    % Verify closed-loop stability with the extracted controller
    try
        closed_loop = feedback(G*K, 1);
        cl_poles = pole(closed_loop);
        
        is_stable = all(real(cl_poles) < 0);
        
        if is_stable
            details = [details, '\nFinal controller stabilizes the plant! Closed-loop poles:\n'];
        else
            details = [details, '\nWARNING: Final controller does not stabilize the plant! Closed-loop poles:\n'];
        end
        
        for i = 1:length(cl_poles)
            if imag(cl_poles(i)) ~= 0
                details = [details, sprintf('  p%d = %.4f + %.4fi\n', i, real(cl_poles(i)), imag(cl_poles(i)))];
            else
                details = [details, sprintf('  p%d = %.4f\n', i, real(cl_poles(i)))];
            end
        end
        
        % If unstable, try final adjustments
        if ~is_stable
            details = [details, 'Attempting final stabilization...\n'];
            
            % Try reducing the gain
            found_stable = false;
            for scale = [0.5, 0.2, 0.1, 0.05, 0.01]
                K_scaled = K * scale;
                closed_loop_scaled = feedback(G * K_scaled, 1);
                
                if all(real(pole(closed_loop_scaled)) < 0)
                    K = K_scaled;
                    details = [details, sprintf('System stabilized by scaling controller gain by %.2f\n', scale)];
                    found_stable = true;
                    break;
                end
            end
            
            if ~found_stable
                details = [details, 'WARNING: Could not stabilize system by gain scaling.\n'];
                
                % Try structure-specific modifications
                details = [details, 'Attempting structure-specific modifications...\n'];
                
                switch structure
                    case 'P'
                        % For P controller, use an extremely conservative gain
                        K = tf(0.01, 1);
                        
                    case 'PI'
                        % For PI, create a controller with minimal integral action
                        K = tf([0.1, 0.001], [1, 0]);
                        
                    case 'PD'
                        % For PD, add more filtering and less derivative action
                        K = tf([0.01, 0.1], [0.1, 1]);
                        
                    case 'PID'
                        % For PID, create a very conservative controller
                        K = tf([0.01, 0.1, 0.001], [0.1, 1, 0]);
                end
                
                % Check if modifications helped
                closed_loop_mod = feedback(G*K, 1);
                is_stable_mod = all(real(pole(closed_loop_mod)) < 0);
                
                if is_stable_mod
                    details = [details, 'Successfully stabilized with structure-specific modifications.\n'];
                else
                    details = [details, 'WARNING: All stabilization attempts failed. Try using pre-stabilization directly.\n'];
                    
                    % Try pre-stabilization as a last resort
                    try
                        [K_prestab, ~] = preStabilize(G, plantInfo);
                        
                        closed_loop_prestab = feedback(G*K_prestab, 1);
                        if all(real(pole(closed_loop_prestab)) < 0)
                            K = K_prestab;
                            details = [details, 'Using pure pre-stabilization controller as final fallback.\n'];
                        else
                            details = [details, 'WARNING: Even pre-stabilization failed. Manual tuning required.\n'];
                            
                            % Use a minimal controller as last resort
                            K = tf(0.001, 1);
                        end
                    catch
                        details = [details, 'Error in pre-stabilization. Manual tuning required.\n'];
                        
                        % Use a minimal controller as last resort
                        K = tf(0.001, 1);
                    end
                end
            end
        end
    catch ME
        details = [details, sprintf('\nError in final verification: %s\n', ME.message)];
    end
    
    % Final stability margin analysis
    try
        [Gm, Pm, Wcg, Wcp] = margin(G*K);
        details = [details, sprintf('\nStability margins:\n')];
        details = [details, sprintf('  Gain margin: %.2f dB at %.4f rad/s\n', 20*log10(Gm), Wcg)];
        details = [details, sprintf('  Phase margin: %.2f deg at %.4f rad/s\n', Pm, Wcp)];
    catch
        details = [details, '\nCould not compute stability margins.\n'];
    end
    
    return;
end
