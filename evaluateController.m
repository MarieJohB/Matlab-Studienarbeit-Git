% EVALUATECONTROLLER Evaluates a controller based on specified criteria
function score = evaluateController(K, G, goal, desired_pm, desired_os, desired_ts, desired_bw)
    % Initialize score
    score = 0;
    
    % Ensure goal is in correct format
    if strcmpi(goal, 'tracking')
        goal = 'Tracking';
    elseif strcmpi(goal, 'disturbance rejection')
        goal = 'Disturbance Rejection';
    elseif strcmpi(goal, 'robustness')
        goal = 'Robustness';
    end
    
    % Compute closed-loop transfer functions
    try
        T = feedback(G*K, 1);      % Complementary sensitivity function
        S = feedback(1, G*K);      % Sensitivity function
        L = G*K;                   % Open-loop transfer function
        
        % Test if closed-loop system is stable
        if any(real(pole(T)) > 0)
            disp('Controller results in unstable closed-loop system.');
            score = -100;  % Unstable system gets a heavily negative score
            return;
        end
    catch ME
        disp(['Error computing closed-loop transfer functions: ', ME.message]);
        score = -50;  % Error indicates likely problems with the controller
        return;
    end
    
    % 1. Stability metrics (30 points maximum)
    try
        % Compute gain and phase margins
        [gm, pm, wgm, wpm] = margin(L);
        
        % Convert gain margin to dB
        gm_dB = 20*log10(gm);
        
        % Phase margin scoring (0-15 points)
        pm_score = 15 * (1 - min(1, abs(pm - desired_pm) / 45));
        
        % Gain margin scoring (0-15 points)
        gm_score = 15 * (1 - min(1, abs(gm_dB - 10) / 10));
        
        % Display stability metrics
        disp(['Phase Margin: ', num2str(pm), '° (desired: ', num2str(desired_pm), '°)']);
        disp(['Gain Margin: ', num2str(gm_dB), ' dB']);
        
        % Add to total score
        score = score + pm_score + gm_score;
    catch ME
        disp(['Warning: Could not compute stability margins: ', ME.message]);
        % Apply penalty for not being able to compute stability margins
        score = score - 15;
    end
    
    % 2. Time-domain performance (40 points maximum)
    try
        % Compute step response
        [y, t] = step(T);
        
        % Extract step response metrics
        info = stepinfo(y, t);
        
        % Calculate key metrics
        actual_overshoot = info.Overshoot;
        actual_settling = info.SettlingTime;
        actual_rise = info.RiseTime;
        
        % Display time domain metrics
        disp(['Overshoot: ', num2str(actual_overshoot), '% (desired: ', num2str(desired_os), '%)']);
        disp(['Settling Time: ', num2str(actual_settling), ' s (desired: ', num2str(desired_ts), ' s)']);
        disp(['Rise Time: ', num2str(actual_rise), ' s']);
        
        % Score overshoot (0-15 points)
        % Better score for being close to desired overshoot
        os_score = 15 * (1 - min(1, abs(actual_overshoot - desired_os) / 50));
        
        % Score settling time (0-15 points)
        % Better score for faster settling
        ts_score = 15 * (1 - min(1, abs(actual_settling - desired_ts) / (2*desired_ts)));
        
        % Score rise time (0-10 points)
        % Better score for faster rise time relative to settling time
        rt_score = 10 * (1 - min(1, actual_rise / desired_ts));
        
        % Add to total score
        score = score + os_score + ts_score + rt_score;
    catch ME
        disp(['Warning: Could not compute time-domain performance metrics: ', ME.message]);
        % Apply penalty for not being able to compute time domain metrics
        score = score - 20;
    end
    
    % 3. Frequency-domain performance (20 points maximum)
    try
        % Calculate bandwidth
        bw = bandwidth(T);
        
        % Display bandwidth
        disp(['Bandwidth: ', num2str(bw), ' rad/s (desired: ', num2str(desired_bw), ' rad/s)']);
        
        % Score bandwidth (0-20 points)
        % Better score for being close to desired bandwidth
        bw_score = 20 * (1 - min(1, abs(bw - desired_bw) / desired_bw));
        
        % Add to total score
        score = score + bw_score;
    catch ME
        disp(['Warning: Could not compute bandwidth: ', ME.message]);
        % Apply penalty for not being able to compute bandwidth
        score = score - 10;
    end
    
    % 4. Goal-specific metrics (10 points maximum)
    try
        switch goal
            case 'Tracking'
                % For tracking, evaluate the reference tracking performance
                % Calculate integral of error for a step reference
                t_sim = linspace(0, 3*desired_ts, 1000);
                r = ones(size(t_sim));
                [y_r, ~] = lsim(T, r, t_sim);
                err_integral = trapz(t_sim, abs(r - y_r));
                
                % Score tracking (0-10 points)
                tracking_score = 10 * exp(-err_integral / 5);
                
                % Add to total score
                score = score + tracking_score;
                
                disp(['Tracking Performance Score: ', num2str(tracking_score), '/10']);
                
            case 'Disturbance Rejection'
                % For disturbance rejection, evaluate the disturbance response
                % Calculate the integral of output for a step disturbance
                t_sim = linspace(0, 3*desired_ts, 1000);
                d = ones(size(t_sim));
                [y_d, ~] = lsim(S, d, t_sim);  % S is sensitivity function
                dist_integral = trapz(t_sim, abs(y_d));
                
                % Score disturbance rejection (0-10 points)
                dist_score = 10 * exp(-dist_integral / 5);
                
                % Add to total score
                score = score + dist_score;
                
                disp(['Disturbance Rejection Score: ', num2str(dist_score), '/10']);
                
            case 'Robustness'
                % For robustness, evaluate the sensitivity peak
                try
                    [Ms, w_ms] = getPeakGain(S);
                    
                    % Good sensitivity peak should be below 2 (6 dB)
                    robust_score = 10 * exp(-(Ms-1) / 1.5);
                    
                    % Add to total score
                    score = score + robust_score;
                    
                    disp(['Sensitivity Peak (Ms): ', num2str(Ms)]);
                    disp(['Robustness Score: ', num2str(robust_score), '/10']);
                catch
                    disp('Could not compute sensitivity peak.');
                    score = score - 5;
                end
        end
    catch ME
        disp(['Warning: Could not compute goal-specific metrics: ', ME.message]);
        % Apply penalty for not being able to compute goal-specific metrics
        score = score - 5;
    end
    
    % Normalize score to be within 0-100 range
    score = max(0, min(100, score));
    
    % Display final score
    disp(['Final Controller Score: ', num2str(score), '/100']);
end