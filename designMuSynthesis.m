function [K, details] = designMuSynthesis(G, structure, options, plantInfo)
% DESIGNMUSYNTHESIS Controller design using robust µ-synthesis approach
%
% Implements an advanced µ-synthesis controller design for systems with uncertainty.
% This is an extension of H-infinity that can handle structured uncertainty more 
% effectively, making it particularly suitable for challenging and unstable plants.
%
% Inputs:
%   G        - Plant transfer function
%   structure - Controller structure ('P', 'PI', 'PD', 'PID')
%   options   - Structure with design parameters:
%     .bandwidth  - Desired bandwidth in rad/s (default: 1)
%     .damping    - Desired damping ratio (default: 0.8)
%     .epsilon    - Filter parameter for D-term (default: 0.1)
%     .robustness - Robustness level ('Low', 'Medium', 'High') (default: 'Medium')
%     .uncertainty - Uncertainty percentage (default: 20)
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
    
    % Define uncertainty based on plant characteristics
    try
        % Convert to state space for modeling uncertainty
        [A, B, C, D] = ssdata(G);
        n = size(A, 1);
        
        % Define uncertainty models based on plant characteristics
        if hasRobustToolbox
            details = [details, 'Creating structured uncertainty model...\n'];
            
            % Create nominal plant model
            G_nom = G;
            
            % Define uncertainty level based on robustness setting
            switch robustness
                case 'Low'
                    delta_gain = uncertaintyPercent / 100;  % Less conservative
                case 'High'
                    delta_gain = 2 * uncertaintyPercent / 100;  % More conservative
                otherwise  % Medium
                    delta_gain = 1.5 * uncertaintyPercent / 100;
            end
            
            % Create input and output uncertainty models
            if plantInfo.isUnstable
                details = [details, 'Unstable plant: Using multiplicative input and output uncertainty models.\n'];
                
                % For unstable plants, use more sophisticated uncertainty modeling
                % Create input and output weights that reflect the uncertainty
                if plantInfo.hasRHPZeros
                    % Non-minimum phase + unstable: more conservative
                    Wout = makeweight(delta_gain/2, omega*2, 2*delta_gain);
                    Win = makeweight(delta_gain/2, omega*2, 2*delta_gain);
                    details = [details, 'Added conservative uncertainty model for non-minimum phase unstable plant.\n'];
                else
                    % Just unstable
                    Wout = makeweight(delta_gain/3, omega*1.5, delta_gain);
                    Win = makeweight(delta_gain/3, omega*1.5, delta_gain);
                end
                
                % Create uncertain plant
                Delta_in = ultidyn('Delta_in', [1 1], 'Bound', 1);
                Delta_out = ultidyn('Delta_out', [1 1], 'Bound', 1);
                
                G_unc = G_nom * (1 + Win*Delta_in) * (1 + Wout*Delta_out);
                
            else
                % For stable plants, simpler uncertainty model is sufficient
                details = [details, 'Stable plant: Using multiplicative input uncertainty model.\n'];
                
                % Create input uncertainty weight
                if plantInfo.hasRHPZeros
                    % Non-minimum phase: more conservative
                    Win = makeweight(delta_gain/3, omega*1.5, delta_gain);
                    details = [details, 'Using conservative uncertainty model for non-minimum phase plant.\n'];
                else
                    Win = makeweight(delta_gain/4, omega, delta_gain/2);
                }
                
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
            
            % Create appropriate frequency-dependent weights
            if plantInfo.isUnstable
                % For unstable plants, more sophisticated weighting
                if plantInfo.hasRHPZeros
                    % Non-minimum phase + unstable: more conservative
                    Wout = tf([delta_gain/2, delta_gain*omega*2], [1, omega*2/delta_gain/2]);
                    Win = tf([delta_gain/2, delta_gain*omega*2], [1, omega*2/delta_gain/2]);
                    details = [details, 'Added conservative uncertainty model for non-minimum phase unstable plant.\n'];
                else
                    % Just unstable
                    Wout = tf([delta_gain/3, delta_gain*omega*1.5], [1, omega*1.5/delta_gain/3]);
                    Win = tf([delta_gain/3, delta_gain*omega*1.5], [1, omega*1.5/delta_gain/3]);
                end
                
                % No need to create actual uncertain plant here as we'll use these weights directly
                % in the generalized plant for D-K iteration
            else
                % For stable plants
                if plantInfo.hasRHPZeros
                    % Non-minimum phase: more conservative
                    Win = tf([delta_gain/3, delta_gain*omega*1.5], [1, omega*1.5/delta_gain/3]);
                    Wout = tf(0.01, 1);  % Minimal output uncertainty for stable plants
                    details = [details, 'Using conservative uncertainty model for non-minimum phase plant.\n'];
                else
                    Win = tf([delta_gain/4, delta_gain*omega/2], [1, omega/delta_gain/4]);
                    Wout = tf(0.01, 1);
                end
            end
            
            G_unc = G;  % Just use nominal plant, weights will be used in generalized plant
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
        
        G_unc = G;  % Use nominal plant
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
            unstable_p = p(real(p) > 0);
            
            if ~isempty(unstable_p)
                fastest_unstable = max(abs(unstable_p));
                Wu = Wu * tf([1, fastest_unstable*2], [1, fastest_unstable/2]);
                details = [details, sprintf('Adjusted control effort weight around unstable pole frequencies (%.3f rad/s).\n', fastest_unstable)];
            end
        end
        
        % Weight on complementary sensitivity function (T = GK/(1+GK))
        % Ensures robustness to unmodeled dynamics and noise rejection
        % Higher robustness requires more roll-off (smaller T) at high frequencies
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
            % Create generalized plant for mu-synthesis
            if plantInfo.isUnstable && plantInfo.hasRHPZeros
                % More complex interconnection for challenging plants
                systemnames = 'G_unc Ws Wu Wt';
                inputvar = '[r; d; n; u]';
                outputvar = '[Ws; Wu; Wt; r-G_unc]';
                input_to_G_unc = '[u]';
                input_to_Ws = '[r-G_unc]';
                input_to_Wu = '[u]';
                input_to_Wt = '[G_unc]';
                
                P = sysic;
            else
                % Standard mixed-sensitivity setup
                systemnames = 'G_unc Ws Wu Wt';
                inputvar = '[r; d; n; u]';
                outputvar = '[Ws; Wu; Wt; r-G_unc]';
                input_to_G_unc = '[u]';
                input_to_Ws = '[r-G_unc]';
                input_to_Wu = '[u]';
                input_to_Wt = '[G_unc]';
                
                P = sysic;
            end
        else
            % Manual creation of generalized plant without Robust Control Toolbox
            details = [details, 'Creating manual generalized plant for D-K iteration.\n'];
            
            % Convert all components to state-space for manual assembly
            [A_G, B_G, C_G, D_G] = ssdata(G);
            [A_Ws, B_Ws, C_Ws, D_Ws] = ssdata(Ws);
            [A_Wu, B_Wu, C_Wu, D_Wu] = ssdata(Wu);
            [A_Wt, B_Wt, C_Wt, D_Wt] = ssdata(Wt);
            
            % Extract dimensions
            n_G = size(A_G, 1);
            n_Ws = size(A_Ws, 1);
            n_Wu = size(A_Wu, 1);
            n_Wt = size(A_Wt, 1);
            
            % Create augmented state-space matrices
            A_P = blkdiag(A_G, A_Ws, A_Wu, A_Wt);
            
            % Input connections [r; d; n; u]
            B_P = [zeros(n_G, 3), B_G;  % G only connects to u
                   B_Ws, zeros(n_Ws, 3);  % Ws connects to error
                   zeros(n_Wu, 3), B_Wu;  % Wu connects to u
                   zeros(n_Wt, 3), zeros(n_Wt, 1)];  % Wt connects to G output
            
            % Modify for correct interconnections
            A_P(n_G+1:n_G+n_Ws, 1:n_G) = B_Ws * -C_G;  % Ws input from error
            A_P(n_G+n_Ws+n_Wu+1:end, 1:n_G) = B_Wt * C_G;  % Wt input from G output
            
            % Output connections [Ws; Wu; Wt; r-G]
            C_P = [zeros(1, n_G), C_Ws, zeros(1, n_Wu), zeros(1, n_Wt);
                   zeros(1, n_G), zeros(1, n_Ws), C_Wu, zeros(1, n_Wt);
                   zeros(1, n_G), zeros(1, n_Ws), zeros(1, n_Wu), C_Wt;
                   -C_G, zeros(1, n_Ws), zeros(1, n_Wu), zeros(1, n_Wt)];
            
            % Feedthrough matrix
            D_P = [0, 0, 0, D_Ws*D_G;
                   0, 0, 0, D_Wu;
                   0, 0, 0, D_Wt*D_G;
                   1, 1, 1, -D_G];
            
            % Create the generalized plant
            P = ss(A_P, B_P, C_P, D_P);
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
        P = augw(G, Ws, Wu, Wt);
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
            
            % Perform mu-synthesis
            [K_mu, CL, mubnds] = musyn(P, 4, 1, opt);
            
            details = [details, sprintf('µ-synthesis completed with %d D-K iterations.\n', length(mubnds))];
            details = [details, sprintf('Final µ-bound: %.4f\n', mubnds(end))];
            
            % Use the resulting controller
            K_robust = K_mu;
        else
            % Implement manual D-K iteration
            details = [details, 'Performing manual D-K iteration approximation...\n'];
            
            % We'll use mixed-sensitivity H-infinity design as approximation
            % to the D-K iteration steps
            
            % Step 1: Initial H-infinity design
            [K1, CL, gamma] = hinfsyn(P, 1, 1);
            details = [details, sprintf('Initial H-infinity design complete with gamma = %.4f\n', gamma)];
            
            % Step 2: Analyze and adjust weights
            % We'll simulate D-K iteration by adjusting weights and redesigning
            
            % Create two more iterations with adjusted weights
            try
                % Extract frequency response of initial controller
                w = logspace(-2, log10(omega*20), 100);
                mag = sigma(CL, w);
                mag = squeeze(mag);
                
                % Find peak frequency
                [peak_mag, idx] = max(mag);
                peak_freq = w(idx);
                
                details = [details, sprintf('Peak magnitude %.4f at frequency %.4f rad/s\n', peak_mag, peak_freq)];
                
                % Adjust weights to focus on critical frequencies
                Ws_adj = Ws * tf([1, peak_freq*0.5], [1, peak_freq*2]);
                Wt_adj = Wt * tf([1, peak_freq*0.5], [1, peak_freq*2]);
                
                % Create adjusted plant
                P_adj = augw(G, Ws_adj, Wu, Wt_adj);
                
                % Second iteration
                [K2, CL2, gamma2] = hinfsyn(P_adj, 1, 1);
                details = [details, sprintf('Second iteration complete with gamma = %.4f\n', gamma2)];
                
                % Third iteration with further adjustment
                if gamma2 < gamma
                    % If improvement, continue in same direction
                    Ws_adj2 = Ws_adj * tf([1, peak_freq*0.2], [1, peak_freq*5]);
                    Wt_adj2 = Wt_adj * tf([1, peak_freq*0.2], [1, peak_freq*5]);
                else
                    % If no improvement, try different approach
                    Ws_adj2 = Ws * tf([1, omega*0.5], [1, omega*5]);
                    Wt_adj2 = Wt * tf([1, omega*2], [1, omega*0.5]);
                }
                
                P_adj2 = augw(G, Ws_adj2, Wu, Wt_adj2);
                
                [K3, CL3, gamma3] = hinfsyn(P_adj2, 1, 1);
                details = [details, sprintf('Third iteration complete with gamma = %.4f\n', gamma3)];
                
                % Select best controller based on gamma
                if gamma3 <= min(gamma, gamma2)
                    K_robust = K3;
                    details = [details, 'Using controller from third iteration.\n'];
                elseif gamma2 <= min(gamma, gamma3)
                    K_robust = K2;
                    details = [details, 'Using controller from second iteration.\n'];
                else
                    K_robust = K1;
                    details = [details, 'Using controller from first iteration.\n'];
                end
            catch ME2
                details = [details, sprintf('Error in iteration: %s\n', ME2.message)];
                details = [details, 'Using controller from initial H-infinity design.\n'];
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
                
                details = [details, sprintf('Reduced controller order from %d to %d\n', 
                          controller_order, keep_states)];
                
                K_robust = K_red;
            else
                details = [details, 'Controller already has low order. No reduction performed.\n'];
            }
        catch ME
            details = [details, sprintf('Order reduction failed: %s\n', ME.message)];
            details = [details, 'Continuing with original controller.\n'];
        end
        
    catch ME
        details = [details, sprintf('Error in controller synthesis: %s\n', ME.message)];
        details = [details, 'Using fallback control design approach...\n'];
        
        % Fallback to standard loop-shaping approach
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
        details = [details, loop_details];
        
        % Skip to the controller extraction step
        if exist('K_robust', 'var')
            % Success, continue with extraction
        else
            % If all else fails, create a very basic controller
            details = [details, 'All synthesis methods failed. Creating basic controller.\n'];
            
            if plantInfo.isUnstable
                % For unstable plants, create a stabilizing controller
                p = plantInfo.poles;
                unstable_poles = p(real(p) > 0);
                
                if ~isempty(unstable_poles)
                    % Simple stabilizing controller: place zeros at unstable poles
                    K_num = 1;
                    K_den = 1;
                    
                    for i = 1:length(unstable_poles)
                        pole_i = unstable_poles(i);
                        
                        if imag(pole_i) ~= 0
                            % Complex pole, need to include conjugate pair
                            if i < length(unstable_poles) && abs(pole_i - conj(unstable_poles(i+1))) < 1e-6
                                % Skip conjugate pair, we'll add both together
                                continue;
                            end
                            
                            % Add complex conjugate pair
                            if imag(pole_i) > 0
                                real_part = real(pole_i);
                                imag_part = imag(pole_i);
                                
                                % Create (s - p)(s - p*)
                                quad_term = [1, -2*real_part, real_part^2 + imag_part^2];
                                
                                % Reflect to LHP: (s + p)(s + p*)
                                stable_quad = [1, 2*real_part, real_part^2 + imag_part^2];
                                
                                K_num = conv(K_num, quad_term);
                                K_den = conv(K_den, stable_quad);
                            end
                        else
                            % Real pole
                            K_num = conv(K_num, [1, -pole_i]);
                            K_den = conv(K_den, [1, abs(pole_i)]);
                        end
                    end
                    
                    % Create the basic controller
                    K_robust = tf(K_num, K_den);
                else
                    % Simple P controller
                    K_robust = tf(1, 1);
                end
            else
                % For stable plants, use a simple controller based on structure
                switch structure
                    case 'P'
                        K_robust = tf(1, 1);
                    case 'PI'
                        K_robust = tf([1, 0.1], [1, 0]);
                    case 'PD'
                        K_robust = tf([0.1, 1], [0.01, 1]);
                    case 'PID'
                        K_robust = tf([0.1, 1, 0.1], [0.01, 1, 0]);
                    otherwise
                        K_robust = tf(1, 1);
                end
            end
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
    
    % Extract controller parameters using frequency response
    try
        w = logspace(-3, log10(omega*10), 200);
        [mag, phase] = bode(K_robust, w);
        mag = squeeze(mag);
        phase = squeeze(phase);
        
        % Create controller based on requested structure
        switch structure
            case 'P'
                % Extract proportional gain (at middle frequencies)
                mid_idx = ceil(length(w)/2);
                Kp = mag(mid_idx);
                
                % Limit gain to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                
                K = tf(Kp, 1);
                details = [details, sprintf('Extracted P controller: Kp = %.4f\n', Kp)];
                
            case 'PI'
                % Extract PI parameters
                
                % Check for integral action (phase approaching -90° at low frequencies)
                has_integrator = (phase(1) < -45);
                
                if has_integrator
                    details = [details, 'Detected integral action in synthesized controller.\n'];
                    
                    % Extract gain at crossover frequency
                    cross_idx = find(mag < 1, 1);
                    if isempty(cross_idx)
                        cross_idx = ceil(length(w)/2);
                    end
                    Kp = mag(cross_idx);
                    
                    % Find corner frequency where integral action starts
                    idx_45 = find(phase > -45, 1);
                    if isempty(idx_45)
                        idx_45 = 2;
                    end
                    w_i = w(idx_45);
                    
                    % Calculate integral gain
                    Ki = Kp * w_i / 5;
                else
                    details = [details, 'No strong integral action detected. Adding appropriate integral term.\n'];
                    
                    % Extract gain at mid frequencies
                    mid_idx = ceil(length(w)/2);
                    Kp = mag(mid_idx);
                    
                    % Add integral action based on bandwidth
                    Ki = Kp * omega / 10;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                        details = [details, 'Reduced integral gain due to plant integrator.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                Ki = min(max(Ki, 0.01), 50);
                
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
                    
                    % Extract parameters
                    Kp = mag(ceil(length(w)/2));  % Mid-frequency gain
                    Kd = Kp / w_d;
                else
                    details = [details, 'No strong derivative action detected. Adding appropriate derivative term.\n'];
                    
                    % Extract gain
                    Kp = mag(ceil(length(w)/2));
                    
                    % Add derivative action based on bandwidth
                    Kd = Kp / omega;
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                        details = [details, 'Reduced derivative action due to RHP zeros.\n'];
                    end
                end
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                Kd = min(max(Kd, 0.01), 50);
                
                % Add filtering for derivative term
                Td = Kd / Kp;
                K = tf([Kd, Kp], [epsilon*Td, 1]);
                details = [details, sprintf('Extracted PD controller: Kp = %.4f, Kd = %.4f\n', Kp, Kd)];
                
            case 'PID'
                % Extract PID parameters
                
                % Check for integral and derivative actions
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
                    details = [details, 'Full PID behavior not detected. Creating appropriate PID controller.\n'];
                    
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
                
                % Limit to reasonable values
                Kp = min(max(Kp, 0.1), 100);
                Ki = min(max(Ki, 0.01), 50);
                Kd = min(max(Kd, 0.01), 50);
                
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
                    
                case 'PI'
                    Ki = K_gain * omega / 10;
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                    end
                    K = tf([K_gain, Ki], [1, 0]);
                    
                case 'PD'
                    Kd = K_gain / omega;
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                    end
                    Td = Kd / K_gain;
                    K = tf([Kd, K_gain], [epsilon*Td, 1]);
                    
                case 'PID'
                    Ki = K_gain * omega / 10;
                    Kd = K_gain / omega;
                    
                    if plantInfo.hasIntegrator
                        Ki = Ki / 2;
                    end
                    
                    if plantInfo.hasRHPZeros
                        Kd = Kd / 2;
                    end
                    
                    K = tf([Kd, K_gain, Ki], [epsilon*Kd, 1, 0]);
                    
                otherwise
                    K = tf(K_gain, 1);
            end
            
            details = [details, 'Created controller using simplistic parameter extraction.\n'];
        catch
            % If everything fails, create a default controller
            switch structure
                case 'P'
                    K = tf(1, 1);
                case 'PI'
                    K = tf([1, 0.1], [1, 0]);
                case 'PD'
                    K = tf([0.1, 1], [epsilon*0.1, 1]);
                case 'PID'
                    K = tf([0.1, 1, 0.1], [epsilon*0.1, 1, 0]);
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
            details = [details, '\nFinal structured controller stabilizes the plant! Closed-loop poles:\n'];
        else
            details = [details, '\nWARNING: Final structured controller does not stabilize the plant! Closed-loop poles:\n'];
        }
        
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
                details = [details, 'WARNING: Could not stabilize system with derived controller.\n'];
                
                % Fallback to a completely different approach
                if plantInfo.isUnstable
                    details = [details, 'Creating backup stabilizing controller for unstable plant.\n'];
                    
                    % Create a simple stabilizing controller directly
                    p = plantInfo.poles;
                    unstable_p = p(real(p) > 0);
                    
                    if ~isempty(unstable_p)
                        % Create stabilizing controller by placing zeros at unstable poles
                        num = 1;
                        den = 1;
                        
                        for i = 1:length(unstable_p)
                            pole_i = unstable_p(i);
                            
                            if imag(pole_i) ~= 0
                                % Skip conjugate pairs, we'll handle them together
                                if i < length(unstable_p) && abs(pole_i - conj(unstable_p(i+1))) < 1e-6
                                    continue;
                                end
                                
                                if imag(pole_i) > 0
                                    real_part = real(pole_i);
                                    imag_part = imag(pole_i);
                                    
                                    % Create zeros at unstable poles
                                    num = conv(num, [1, -2*real_part, real_part^2 + imag_part^2]);
                                    
                                    % Create stable poles
                                    den = conv(den, [1, 2*abs(real_part), abs(real_part)^2 + imag_part^2]);
                                end
                            else
                                % Real pole
                                num = conv(num, [1, -pole_i]);
                                den = conv(den, [1, 2*abs(pole_i)]);
                            end
                        end
                        
                        % Create stabilizing controller
                        K_stab = tf(0.1*num, den);
                        
                        % Check stability
                        closed_loop_stab = feedback(G*K_stab, 1);
                        if all(real(pole(closed_loop_stab)) < 0)
                            K = K_stab;
                            details = [details, 'Basic stabilizing controller created successfully.\n'];
                            
                            % Convert to requested structure if needed
                            if strcmpi(structure, 'PID')
                                % Create an approximate PID
                                try
                                    [num_stab, den_stab] = tfdata(K_stab, 'v');
                                    Kp = 0.1;
                                    Ki = 0.01;
                                    Kd = 0.5;
                                    K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                                    
                                    % Make sure it's still stabilizing
                                    if ~all(real(pole(feedback(G*K, 1))) < 0)
                                        K = K_stab;  % Revert if not stable
                                    end
                                catch
                                    K = K_stab;  // Keep original if conversion fails
                                end
                            end
                        else
                            details = [details, 'Even basic stabilizing controller failed. This system is exceptionally challenging.\n'];
                        end
                    end
                else
                    // For stable plants, use a conservative controller for safety
                    switch structure
                        case 'P'
                            K = tf(0.1, 1);
                        case 'PI'
                            K = tf([0.1, 0.01], [1, 0]);
                        case 'PD'
                            K = tf([0.05, 0.1], [0.01, 1]);
                        case 'PID'
                            K = tf([0.05, 0.1, 0.01], [0.01, 1, 0]);
                        otherwise
                            K = tf(0.1, 1);
                    end
                    
                    details = [details, 'Using very conservative controller parameters.\n'];
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

% Helper function to create a formatted string with plant information
function infoStr = getPlantInfoString(plantInfo)
    % Initialize output string
    infoStr = '';
    
    % Add stability information
    if plantInfo.isUnstable
        infoStr = [infoStr, 'Unstable, '];
    else
        infoStr = [infoStr, 'Stable, '];
    end
    
    % Add integrator information
    if plantInfo.hasIntegrator
        infoStr = [infoStr, 'Has integrator, '];
    end
    
    % Add RHP zeros information
    if plantInfo.hasRHPZeros
        infoStr = [infoStr, 'Non-minimum phase, '];
    end
    
    % Add delay information
    if plantInfo.hasDelay
        infoStr = [infoStr, 'Has delay, '];
    end
    
    % Add order information
    if plantInfo.isHighOrder
        infoStr = [infoStr, 'High-order, '];
    else
        infoStr = [infoStr, 'Low-order, '];
    end
    
    % Add DC gain if available
    if ~isnan(plantInfo.dcGain) && ~isinf(plantInfo.dcGain)
        infoStr = [infoStr, sprintf('DC gain=%.3g', plantInfo.dcGain)];
    else
        infoStr = [infoStr, 'Infinite DC gain'];
    end
end