function K = get_user_controller(G)
% get_user_controller
% ---------------------
% This function creates a MATLAB app for selecting a controller type and
% entering its parameters with a live HTML preview. The input dialog layout,
% including the preview area, labels, input fields and buttons,
% is exactly modeled after the get_user_transfer_function function.
%
% New: Now includes an "Auto-Tuning" button that creates an optimized controller
% based on various design methods.

    % Initialize K to empty to ensure it always has a value
    K = [];
    
    % Controller selection menu using a uifigure with buttons
    mainFig = uifigure('Name', 'Controller Selection', 'Position', [100 100 600 400]);
    uilabel(mainFig, 'Text', 'Select a Controller:', ...
        'Position', [200 350 200 30], 'HorizontalAlignment', 'center', 'FontSize', 16);
    
    % Create buttons for each controller type
    controllers = {'P','PI','PD','PID','PT1','PIT1','I2','PIDT1','Custom'};
    btnPositions = [50 250; 170 250; 290 250; 410 250; 50 170; 170 170; 290 170; 410 170; 250 90];
    for i = 1:numel(controllers)
        uibutton(mainFig, 'push', 'Text', controllers{i}, ...
            'Position', [btnPositions(i,1) btnPositions(i,2) 100 40], ...
            'ButtonPushedFcn', @(~,~) openInputDialog(controllers{i}));
    end
    
    % Add the Auto-Tuning button
    autoTuneBtn = uibutton(mainFig, 'push', 'Text', 'Auto-Tuning', ...
        'Position', [185 30 230 50], 'FontSize', 14, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.3 0.6 0.9], 'FontColor', 'white', ...
        'ButtonPushedFcn', @(~,~) openAutoTuningDialog());
    
    % Handle figure close request
    mainFig.CloseRequestFcn = @handleCloseRequest;
    
    % Block execution until a valid controller is confirmed
    uiwait(mainFig);
    
    %% Nested Functions
    
    % Handle close request for main figure
    function handleCloseRequest(~,~)
        % Ensure K has a value (already initialized to [])
        % This will be the return value when user closes without selecting
        uiresume(mainFig);
        delete(mainFig);
    end
    
    % Opens the Auto-Tuning dialog
    function openAutoTuningDialog()
        % This function creates an improved auto-tuning dialog for controller design
        % It provides a more organized UI with clear sections, fixed alignment, tooltips,
        % and extended methods including H-infinity and LQG

        % Use the G plant model passed as parameter to the outer function
        % We're simply checking if G is empty - not checking nargin which would be wrong in this context
        if isempty(G)
             uialert(mainFig, 'Please define plant model G(s) first.', 'Plant Missing');
             return;
        end
        
        % Check stability and order of G
        isStable = all(real(pole(G)) < 0);
        plantOrder = length(pole(G));
        
        % Create a figure with improved styling
        autoTuneFig = uifigure('Name', 'Automatic Controller Design', 'Position', [300 150 700 700]);
        autoTuneFig.Color = [0.95 0.95 0.97]; % Light gray-blue background
        
        % Add a professional-looking title and plant info
        titlePanel = uipanel(autoTuneFig, 'Position', [10 630 680 60], 'BackgroundColor', [0.2 0.4 0.7]);
        titleLabel = uilabel(titlePanel, 'Text', 'Controller Auto-Tuning', ...
            'Position', [20 15 640 30], 'FontSize', 20, 'FontWeight', 'bold', ...
            'FontColor', 'white', 'HorizontalAlignment', 'center');
        
        % Plant info display with improved styling
        if isStable
            stabilityText = 'stable';
            stabilityColor = [0.2 0.7 0.3]; % Green for stable
        else
            stabilityText = 'unstable';
            stabilityColor = [0.7 0.3 0.2]; % Red for unstable
        end
        
        plantInfoText = sprintf('Plant G(s): %s, Order %d', stabilityText, plantOrder);
        plantInfoLabel = uilabel(autoTuneFig, 'Text', plantInfoText, ...
            'Position', [20 600 660 20], 'FontWeight', 'bold', ...
            'FontColor', stabilityColor);
        
        % Create main panels with improved styling
        designPanel = uipanel(autoTuneFig, 'Title', 'Design Methods', ...
            'Position', [10 420 680 170], 'TitlePosition', 'centertop', ...
            'FontWeight', 'bold', 'FontSize', 14);
        
        controllerPanel = uipanel(autoTuneFig, 'Title', 'Controller Configuration', ...
            'Position', [10 330 680 80], 'TitlePosition', 'centertop', ...
            'FontWeight', 'bold', 'FontSize', 14);
        
        performancePanel = uipanel(autoTuneFig, 'Title', 'Performance Requirements', ...
            'Position', [10 190 680 130], 'TitlePosition', 'centertop', ...
            'FontWeight', 'bold', 'FontSize', 14);
        
        resultsPanel = uipanel(autoTuneFig, 'Title', 'Results', ...
            'Position', [10 10 680 170], 'TitlePosition', 'centertop', ...
            'FontWeight', 'bold', 'FontSize', 14);
        
        % Methods selection with better layout (2 columns, no duplicates)
        methods = {
            'Ziegler-Nichols (Oscillation)', isApplicable(G, 'ZN-Oscillation'),
            'Ziegler-Nichols (Step)', isApplicable(G, 'ZN-Step'),
            'Aström', isApplicable(G, 'Astrom'),
            'CHR (with 0% Overshoot)', isApplicable(G, 'CHR'),
            'Cohen-Coon', isApplicable(G, 'Cohen-Coon'),
            'Loop-Shaping', isApplicable(G, 'Loop-Shaping'),
            'IMC (Internal Model Control)', isApplicable(G, 'IMC'),
            'MIGO (M-constrained Integral Gain Optimization)', isApplicable(G, 'MIGO'),
            'H-infinity', isApplicable(G, 'H-infinity'),
            'LQG (Linear-Quadratic-Gaussian)', isApplicable(G, 'LQG')
        };
        
        % Calculate rows and columns for methods layout
        numMethods = size(methods, 1);
        numRows = ceil(numMethods / 2);
        
        % Create checkboxes in a grid layout
        methodCheckboxes = [];
        for i = 1:numMethods
            row = ceil(i / 2);
            col = mod(i-1, 2);
            
            % Calculate position
            x = 20 + col * 340;
            y = 140 - row * 25;
            
            % Create checkbox with improved styling
            cb = uicheckbox(designPanel, 'Text', methods{i, 1}, ...
                'Position', [x y 320 22], 'Value', methods{i, 2}, ...
                'FontSize', 12);
            
            % Gray out and disable if not applicable
            if ~methods{i, 2}
                cb.Enable = 'off';
                cb.Value = false;
                cb.Text = [methods{i, 1}, ' (not applicable)'];
                cb.FontColor = [0.5 0.5 0.5];
            end
            
            methodCheckboxes = [methodCheckboxes; cb];
        end
        
        % Add a "What are these methods?" button with information tooltip
        infoBtn = uibutton(designPanel, 'Text', 'Method Information', ...
            'Position', [550 20 110 30], 'ButtonPushedFcn', @showMethodInfo);
        
        % Controller structure dropdown with improved styling
        structureLabel = uilabel(controllerPanel, 'Text', 'Controller Structure:', ...
            'Position', [20 30 150 22], 'FontSize', 12);
        structureDropdown = uidropdown(controllerPanel, ...
            'Items', {'P', 'PI', 'PD', 'PID'}, ...
            'Position', [170 30 150 22], 'Value', 'PID', ...
            'FontSize', 12);
        
        % D-Term filter parameter with improved styling
        filterLabel = uilabel(controllerPanel, 'Text', 'D-Term Filter (epsilon):', ...
            'Position', [350 30 170 22], 'FontSize', 12);
        filterField = uieditfield(controllerPanel, 'numeric', ...
            'Position', [520 30 100 22], 'Value', 0.1, ...
            'Limits', [0.001 0.5], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
            'FontSize', 12);
        
        % Performance requirements with improved layout
        % Left column
        phaseMarginLabel = uilabel(performancePanel, 'Text', 'Phase Margin (deg):', ...
            'Position', [20 80 150 22], 'FontSize', 12);
        phaseMarginField = uieditfield(performancePanel, 'numeric', ...
            'Position', [170 80 100 22], 'Value', 45, ...
            'Limits', [10 80], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
            'FontSize', 12);
        
        gainMarginLabel = uilabel(performancePanel, 'Text', 'Gain Margin (dB):', ...
            'Position', [20 50 150 22], 'FontSize', 12);
        gainMarginField = uieditfield(performancePanel, 'numeric', ...
            'Position', [170 50 100 22], 'Value', 8, ...
            'Limits', [3 20], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
            'FontSize', 12);
        
        bandwidthLabel = uilabel(performancePanel, 'Text', 'Bandwidth (rad/s):', ...
            'Position', [20 20 150 22], 'FontSize', 12);
        bandwidthField = uieditfield(performancePanel, 'numeric', ...
            'Position', [170 20 100 22], 'Value', 1, ...
            'Limits', [0.01 100], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
            'FontSize', 12);
        
        % Right column
        robustnessLabel = uilabel(performancePanel, 'Text', 'Robustness:', ...
            'Position', [350 80 100 22], 'FontSize', 12);
        robustnessDropdown = uidropdown(performancePanel, ...
            'Items', {'Low', 'Medium', 'High'}, ...
            'Position', [450 80 100 22], 'Value', 'Medium', ...
            'FontSize', 12);
        
        overshootLabel = uilabel(performancePanel, 'Text', 'Overshoot (%):', ...
            'Position', [350 50 100 22], 'FontSize', 12);
        overshootField = uieditfield(performancePanel, 'numeric', ...
            'Position', [450 50 100 22], 'Value', 10, ...
            'Limits', [0 50], 'LowerLimitInclusive', true, 'UpperLimitInclusive', true, ...
            'FontSize', 12);
        
        goalLabel = uilabel(performancePanel, 'Text', 'Optimization Goal:', ...
            'Position', [350 20 120 22], 'FontSize', 12);
        goalDropdown = uidropdown(performancePanel, ...
            'Items', {'Tracking', 'Disturbance Rejection', 'Robustness'}, ...
            'Position', [450 20 170 22], 'Value', 'Tracking', ...
            'FontSize', 12);
        
        % Add a "What's a good score?" button with information
        scoreInfoBtn = uibutton(performancePanel, 'Text', 'Score Info', ...
            'Position', [580 20 80 22], 'ButtonPushedFcn', @showScoreInfo);
        
        % Results area with improved styling
        resultArea = uitextarea(resultsPanel, ...
            'Position', [20 60 640 80], ...
            'Value', 'Auto-tuning results will be displayed here.', ...
            'Editable', 'off', 'FontSize', 12);
        
        % Buttons with improved styling and layout
        startBtn = uibutton(resultsPanel, 'push', 'Text', 'Start Auto-Tuning', ...
            'Position', [120 20 160 30], 'FontSize', 14, 'FontWeight', 'bold', ...
            'BackgroundColor', [0.3 0.8 0.3], 'FontColor', 'white', ...
            'ButtonPushedFcn', @startAutoTuning);
        
        applyBtn = uibutton(resultsPanel, 'push', 'Text', 'Apply', ...
            'Position', [300 20 120 30], 'FontSize', 14, ...
            'Enable', 'off', 'BackgroundColor', [0.3 0.5 0.8], 'FontColor', 'white', ...
            'ButtonPushedFcn', @confirmController);
        
        cancelBtn = uibutton(resultsPanel, 'push', 'Text', 'Cancel', ...
            'Position', [440 20 120 30], 'FontSize', 14, ...
            'BackgroundColor', [0.8 0.3 0.3], 'FontColor', 'white', ...
            'ButtonPushedFcn', @cancelDialog);
        
        % Handle figure close request
        autoTuneFig.CloseRequestFcn = @cancelDialog;
        
        % Storage for the optimized controller
        bestController = [];
        controllerDetails = '';
        
        % ========================
        % Function Implementations
        % ========================
        
        % Show method information
        function showMethodInfo(~, ~)
            methodInfoFig = uifigure('Name', 'Controller Design Methods', 'Position', [350 250 600 500]);
            
            methodInfo = uitextarea(methodInfoFig, 'Position', [20 60 560 420], 'Editable', 'off');
            methodInfo.Value = {
                'Controller Design Methods:', 
                '-------------------------',
                '1. Ziegler-Nichols (Oscillation): Uses critical gain and frequency to determine controller parameters.',
                '   • Best for: Systems that can be safely brought to the stability boundary for testing.',
                '   • Limitations: Can result in aggressive control with high overshoot.',
                '',
                '2. Ziegler-Nichols (Step): Uses open-loop step response characteristics to determine parameters.',
                '   • Best for: Stable plants with S-shaped step response.',
                '   • Limitations: Not suitable for integrating or unstable plants.',
                '',
                '3. Aström: Improved version of Ziegler-Nichols with better robustness.',
                '   • Best for: Similar to Ziegler-Nichols, but when less overshoot is desired.',
                '',
                '4. CHR (Chien-Hrones-Reswick): Designed for either servo or regulator problems with focus on overshoot.',
                '   • Best for: When specific overshoot constraints exist.',
                '',
                '5. Cohen-Coon: Designed specifically for plants with significant dead time.',
                '   • Best for: Processes with large transport delays.',
                '',
                '6. Loop-Shaping: Frequency-domain method focused on achieving desired loop shape.',
                '   • Best for: When specific bandwidth and phase margin are required.',
                '',
                '7. IMC (Internal Model Control): Based on process model and desired closed-loop response.',
                '   • Best for: Good disturbance rejection with specified settling time.',
                '',
                '8. MIGO: Maximizes integral gain while satisfying robustness constraints.',
                '   • Best for: Good balance between performance and robustness.',
                '',
                '9. H-infinity: Robust control method based on minimizing worst-case scenarios.',
                '   • Best for: Systems with significant uncertainties or disturbances.',
                '',
                '10. LQG: Combines linear quadratic regulator with Kalman filtering.',
                '    • Best for: Optimal control in presence of Gaussian noise.'
            };
            
            closeBtn = uibutton(methodInfoFig, 'push', 'Text', 'Close', ...
                'Position', [250 20 100 30], 'FontSize', 12, ...
                'ButtonPushedFcn', @(~,~) close(methodInfoFig));
        end
        
        % Show score information
        function showScoreInfo(~, ~)
            scoreInfoFig = uifigure('Name', 'Controller Scoring System', 'Position', [350 300 500 400]);
            
            scoreInfo = uitextarea(scoreInfoFig, 'Position', [20 60 460 320], 'Editable', 'off');
            scoreInfo.Value = {
                'Controller Scoring System:', 
                '-------------------------',
                'The auto-tuning process evaluates controllers using a point system (0-100):',
                '',
                '• Excellent: 80+ points',
                '  - Meets all design requirements with excellent balance',
                '  - Optimal time and frequency domain performance',
                '  - Strong robustness to disturbances and model uncertainties',
                '',
                '• Good: 60-80 points',
                '  - Meets all key requirements with good compromise',
                '  - Stable with good time response characteristics',
                '  - Suitable phase and gain margins (>45° and >6dB)',
                '',
                '• Acceptable: 40-60 points',
                '  - Meets minimum requirements but with compromises',
                '  - May have some overshoot or longer settling time',
                '  - Sufficient but not optimal robustness',
                '',
                '• Poor: 20-40 points',
                '  - Meets basic stability requirements but performance is lacking',
                '  - May have excessive overshoot or very slow response',
                '  - Limited robustness to disturbances or model variations',
                '',
                '• Unacceptable: <20 points',
                '  - May be stable but with poor performance metrics',
                '  - Fails to meet several key design requirements',
                '  - Not recommended for implementation',
                '',
                'Note: A negative score indicates an unstable controller that should not be used.'
            };
            
            closeBtn = uibutton(scoreInfoFig, 'push', 'Text', 'Close', ...
                'Position', [200 20 100 30], 'FontSize', 12, ...
                'ButtonPushedFcn', @(~,~) close(scoreInfoFig));
        end
        
        % Function to check if a design method is applicable to the plant
        function applicable = isApplicable(G, methodName)
            % Get plant information
            [z, p, k] = zpkdata(G, 'v');
            isStable = all(real(p) < 0);
            hasIntegrator = any(abs(p) < 1e-6);
            hasRHPZero = any(real(z) > 0);
            plantOrder = length(p);
            
            % Check applicability based on plant characteristics
            switch methodName
                case 'ZN-Oscillation'
                    % Works for most plants that can be stabilized
                    applicable = true;
                case 'ZN-Step'
                    % Only for stable plants without integrators
                    applicable = isStable && ~hasIntegrator;
                case 'Astrom'
                    % For stable plants
                    applicable = isStable;
                case 'CHR'
                    % For stable plants
                    applicable = isStable;
                case 'Cohen-Coon'
                    % For stable plants
                    applicable = isStable;
                case 'Loop-Shaping'
                    % All plants
                    applicable = true;
                case 'IMC'
                    % Best for stable, minimum-phase plants but can be adapted
                    applicable = true;
                case 'MIGO'
                    % For stable plants primarily
                    applicable = true;
                case 'H-infinity'
                    % All plants, especially good with uncertainties
                    applicable = true;
                case 'LQG'
                    % All plants, requires state-space representation
                    applicable = true;
                otherwise
                    applicable = false;
            end
        end
        
        % Start Auto-Tuning
        function startAutoTuning(~, ~)
            % Disable the Start button during calculation
            startBtn.Enable = 'off';
            startBtn.Text = 'Processing...';
            drawnow;
            
            % Parameters structure
            options = struct();
            options.epsilon = filterField.Value;
            options.phaseMargin = phaseMarginField.Value;
            options.gainMargin = gainMarginField.Value;
            options.bandwidth = bandwidthField.Value;
            options.settlingTime = 5;  % Default value
            options.robustness = robustnessDropdown.Value;
            options.overshoot = overshootField.Value;
            options.goal = goalDropdown.Value;
            
            % Controller structure
            structure = structureDropdown.Value;
            
            % Collect selected methods
            selectedMethods = {};
            for i = 1:length(methodCheckboxes)
                if methodCheckboxes(i).Value && methodCheckboxes(i).Enable == 'on'
                    methodName = methodCheckboxes(i).Text;
                    % Remove "(not applicable)" if present
                    methodName = strrep(methodName, ' (not applicable)', '');
                    selectedMethods{end+1} = methodName;
                end
            end
            
            if isempty(selectedMethods)
                uialert(autoTuneFig, 'Please select at least one design method.', 'No Method Selected');
                startBtn.Enable = 'on';
                startBtn.Text = 'Start Auto-Tuning';
                return;
            end
            
            % Update status
            resultArea.Value = 'Auto-tuning in progress... This may take a moment.';
            drawnow;
            
            % Storage for best controller
            bestScore = -Inf;
            bestController = [];
            controllerDetails = '';
            
            try
                % Process each method
                for i = 1:length(selectedMethods)
                    methodName = selectedMethods{i};
                    resultArea.Value = ['Testing method: ' methodName '...'];
                    drawnow;
                    
                    try
                        % Make sure design_controller_auto function exists and is working
                        % This is a check for the function, remove if you're sure it exists
                        if ~exist('design_controller_auto', 'file')
                            error('Function design_controller_auto not found. Please make sure it is in your MATLAB path.');
                        end
                        
                        % Design controller using selected method
                        [K_candidate, details] = design_controller_auto(G, methodName, structure, options);
                        
                        % Extract score from details
                        score = extractScoreFromDetails(details);
                        
                        % Update the result area with current method results
                        resultArea.Value{end+1} = sprintf('Method: %s, Score: %.2f', methodName, score);
                        drawnow;
                        
                        % Check if this is the best controller so far
                        if score > bestScore
                            bestScore = score;
                            bestController = K_candidate;
                            scoreLabel = getScoreLabel(score);
                            controllerDetails = sprintf('Method: %s\nScore: %.2f (%s)\n\n%s', methodName, score, scoreLabel, details);
                        end
                    catch ME
                        % Log error for this method with detailed stack trace
                        disp(['Error with method ' methodName ': ' ME.message]);
                        disp(getReport(ME, 'extended'));
                        resultArea.Value{end+1} = sprintf('Error with %s: %s', methodName, ME.message);
                        drawnow;
                    end
                end
                
                % Display final results
                if ~isempty(bestController)
                    scoreLabel = getScoreLabel(bestScore);
                    resultArea.Value = {
                        'Auto-tuning completed!',
                        sprintf('Best controller designed using: %s', extractMethodName(controllerDetails)),
                        sprintf('Score: %.2f (%s)', bestScore, scoreLabel),
                        '',
                        'Check console for detailed results or click "Apply" to use this controller.'
                    };
                    
                    % Enable Apply button
                    applyBtn.Enable = 'on';
                    
                    % Display the controller in the command window
                    disp('====== Auto-Tuning Results ======');
                    disp(controllerDetails);
                    
                    % Get controller details for display
                    [num, den] = tfdata(bestController, 'v');
                    resultArea.Value{end+1} = '';
                    resultArea.Value{end+1} = 'Controller Transfer Function:';
                    resultArea.Value{end+1} = sprintf('Numerator: [%s]', num2str(num, '%.4g '));
                    resultArea.Value{end+1} = sprintf('Denominator: [%s]', num2str(den, '%.4g '));
                else
                    resultArea.Value = 'Auto-tuning failed. Please try different methods or parameters.';
                end
            catch ME
                % Handle general errors with detailed stack trace
                disp(['General error during auto-tuning: ' ME.message]);
                disp(getReport(ME, 'extended'));
                
                resultArea.Value = {
                    'Error during auto-tuning:',
                    ME.message,
                    '',
                    'Please try different parameters or methods.'
                };
            end
            
            % Re-enable the Start button
            startBtn.Enable = 'on';
            startBtn.Text = 'Start Auto-Tuning';
        end
        
        % Extract score from controller details
        function score = extractScoreFromDetails(details)
            % Default score if not found
            score = 0;
            
            % Look for the score information in the details string
            try
                matches = regexp(details, 'Controller Score: (\d+\.\d+)/100', 'tokens');
                if ~isempty(matches) && ~isempty(matches{1})
                    score = str2double(matches{1}{1});
                end
            catch
                % If there's any error, return default score
                score = 0;
            end
        end
        
        % Helper function to extract method name from details
        function methodName = extractMethodName(details)
            lines = strsplit(details, '\n');
            if ~isempty(lines) && length(lines) > 0
                firstLine = lines{1};
                parts = strsplit(firstLine, ':');
                if length(parts) > 1
                    methodName = strtrim(parts{2});
                else
                    methodName = 'Unknown';
                end
            else
                methodName = 'Unknown';
            end
        end
        
        % Helper function to get a textual label for a score
        function label = getScoreLabel(score)
            if score < 0
                label = 'Unstable';
            elseif score < 20
                label = 'Unacceptable';
            elseif score < 40
                label = 'Poor';
            elseif score < 60
                label = 'Acceptable';
            elseif score < 80
                label = 'Good';
            else
                label = 'Excellent';
            end
        end
        
        % Confirm Controller Function
        function confirmController(~, ~)
            if ~isempty(bestController)
                % Set K to the best controller and close dialog
                K = bestController;
                uiresume(mainFig);
                close(autoTuneFig);
                close(mainFig);
            else
                uialert(autoTuneFig, 'No valid controller available. Please run auto-tuning first.', 'No Controller');
            end
        end
        
        % Cancel Function
        function cancelDialog(~, ~)
            % Close auto-tuning window without selecting a controller
            close(autoTuneFig);
        end
    end

    % Modal dialog for controller parameter input
    function openInputDialog(ctrlType)
        % Determine field labels based on controller type
        switch ctrlType
            case 'P'
                fieldLabels = {'Enter Kp:'};
            case 'PI'
                fieldLabels = {'Enter Kp:', 'Enter Ki:'};
            case 'PD'
                fieldLabels = {'Enter Kp:', 'Enter Kd:'};
            case 'PID'
                fieldLabels = {'Enter Kp:', 'Enter Ki:', 'Enter Kd:'};
            case 'PT1'
                fieldLabels = {'Enter Kp:', 'Enter Time Constant T:'};
            case 'PIT1'
                fieldLabels = {'Enter Kp:', 'Enter Ki:', 'Enter Time Constant T:'};
            case 'I2'
                fieldLabels = {'Enter Ki:'};
            case 'PIDT1'
                fieldLabels = {'Enter Kp:', 'Enter Ki:', 'Enter Kd:', 'Enter Time Constant T:'};
            case 'Custom'
                fieldLabels = {'Enter numerator coefficients [b0 b1 ...]:', ...
                               'Enter denominator coefficients [a0 a1 ...]:'};
        end
        numFields = numel(fieldLabels);
        
        % Adjust figure height based on number of input fields (mimic original layout)
        if numFields <= 2
            figHeight = 350;
        else
            figHeight = 350 + (numFields - 2)*50;
        end
        
        % Create a modal UI figure similar to get_user_transfer_function
        dlg = uifigure('Name', [ctrlType ' Controller Input'], ...
                       'Position', [100 100 500 figHeight], 'WindowStyle', 'modal');
                   
        % Create UIHTML component for the preview at the top of the window.
        % Initially, display the symbolic preview.
        previewHTML = uihtml(dlg, 'Position', [20 figHeight-70 460 50], ...
            'HTMLSource', getSymbolicPreviewHTML(ctrlType));
        
        % Create labels and input fields for controller parameters
        fields = cell(1, numFields);
        for j = 1:numFields
            yLabel = figHeight - (130 + (j-1)*60);
            yEdit  = figHeight - (160 + (j-1)*60);
            uilabel(dlg, 'Position', [20 yLabel 250 22], 'Text', fieldLabels{j});
            fields{j} = uieditfield(dlg, 'text', 'Position', [20 yEdit 460 22]);
        end
        
        % Create three centered buttons at the bottom (same layout as transfer function input)
        uibutton(dlg, 'push', 'Text', 'Preview', ...
            'Position', [80 60 100 30], 'ButtonPushedFcn', @(~,~) previewCallbackK());
        uibutton(dlg, 'push', 'Text', 'Confirm', ...
            'Position', [200 60 100 30], 'ButtonPushedFcn', @(~,~) confirmCallbackK());
        uibutton(dlg, 'push', 'Text', 'Cancel', ...
            'Position', [320 60 100 30], 'ButtonPushedFcn', @(~,~) cancelCallbackK());
            
        % Handle figure close request
        dlg.CloseRequestFcn = @cancelCallbackK;
        
        % Callback for Preview: update the HTML preview area
        function previewCallbackK()
            % If any input field is empty, show the symbolic preview
            if any(cellfun(@(x) isempty(x.Value), fields))
                previewHTML.HTMLSource = getSymbolicPreviewHTML(ctrlType);
                return;
            end
            % Otherwise, try to compute numeric transfer function preview
            try
                Ktemp = computeController(ctrlType, fields);
                % Extract numerator and denominator from the transfer function (assumed SISO)
                numCell = Ktemp.Numerator;
                denCell = Ktemp.Denominator;
                if isempty(numCell) || isempty(denCell)
                    previewHTML.HTMLSource = '<html><body style="font-size:16px; text-align:center;">Preview: Invalid input.</body></html>';
                    return;
                end
                numVec = numCell{1};
                denVec = denCell{1};
                % Convert polynomial coefficients to HTML-formatted strings without asterisks
                numHTML = polyToHTMLString(numVec);
                denHTML = polyToHTMLString(denVec);
                % Create HTML fraction with "Preview:" label and centered fraction
                formulaHTML = [...
                    '<html><head><style>',...
                        '.container { display: inline-flex; align-items: center; justify-content: center; width: 100%; }',...
                        '.label { font-size:16px; margin-right:10px; }',...
                        '.fraction { display:inline-block; text-align:center; }',...
                        '.fraction .num { display:block; border-bottom:1px solid black; padding-bottom:2px; }',...
                        '.fraction .den { display:block; padding-top:2px; }',...
                    '</style></head><body>',...
                        '<div class="container">',...
                            '<div class="label">Preview:</div>',...
                            '<div class="fraction">',...
                                '<span class="num">', numHTML, '</span>',...
                                '<span class="den">', denHTML, '</span>',...
                            '</div>',...
                        '</div>',...
                    '</body></html>'];
                previewHTML.HTMLSource = formulaHTML;
            catch ME
                disp(['Preview error: ' ME.message]);
                previewHTML.HTMLSource = '<html><body style="font-size:16px; text-align:center;">Preview: Invalid input.</body></html>';
            end
        end
        
        % Callback for Confirm: validate input, create transfer function and close dialog & main figure
        function confirmCallbackK()
            try
                K = computeController(ctrlType, fields);
                % Resume and close safely
                uiresume(dlg);
                delete(dlg);
                uiresume(mainFig);
                delete(mainFig);
            catch ME
                disp(['Error confirming controller: ' ME.message]);
                uialert(dlg, 'Invalid input. Please enter valid parameter values.', 'Error');
            end
        end
        
        % Callback for Cancel: close the dialog and return to the controller selection window
        function cancelCallbackK()
            uiresume(dlg);
            delete(dlg);
        end
        
        uiwait(dlg);
    end

    % Helper function: Computes the controller transfer function based on controller type and inputs
    function Ktemp = computeController(ctrlType, fields)
        switch ctrlType
            case 'P'
                Kp = convertValue(fields{1}.Value);
                Ktemp = tf(Kp, 1);
            case 'PI'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                Ktemp = tf([Kp, Ki], [1, 0]);
            case 'PD'
                Kp = convertValue(fields{1}.Value);
                Kd = convertValue(fields{2}.Value);
                Ktemp = tf([Kd, Kp], 1);
            case 'PID'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                Kd = convertValue(fields{3}.Value);
                Ktemp = tf([Kd, Kp, Ki], [1, 0]);
            case 'PT1'
                Kp = convertValue(fields{1}.Value);
                T  = convertValue(fields{2}.Value);
                Ktemp = tf(Kp, [T, 1]);
            case 'PIT1'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                T  = convertValue(fields{3}.Value);
                Ktemp = tf([Kp*Ki, Kp], [T, 1, 0]);
            case 'I2'
                Ki = convertValue(fields{1}.Value);
                Ktemp = tf(Ki, [1, 0, 0]);
            case 'PIDT1'
                Kp = convertValue(fields{1}.Value);
                Ki = convertValue(fields{2}.Value);
                Kd = convertValue(fields{3}.Value);
                T  = convertValue(fields{4}.Value);
                Ktemp = tf([Kd, Kp, Ki], [T, 1, 0]);
            case 'Custom'
                numStr = strrep(fields{1}.Value, ',', '.');
                denStr = strrep(fields{2}.Value, ',', '.');
                numVec = str2num(numStr);  %#ok<ST2NM>
                denVec = str2num(denStr);  %#ok<ST2NM>
                if isempty(numVec) || isempty(denVec) || any(isnan(numVec)) || any(isnan(denVec))
                    error('Invalid input');
                end
                Ktemp = tf(numVec, denVec);
            otherwise
                error('Unknown controller type.');
        end
    end

    % Helper function: Converts an input string (allowing decimal commas) to a number.
    % If the input string contains spaces, it takes only the first token.
    function num = convertValue(str)
        % Replace commas with periods
        str = strrep(str, ',', '.');
        % Tokenize the string at whitespace and take only the first token
        token = strtok(str);
        num = str2double(token);
        if isnan(num)
            error('Conversion error');
        end
    end

    % Helper function: Converts a polynomial coefficient vector to an HTML-formatted
    % polynomial string without an asterisk between the coefficient and 's', using <sup> tags for exponents.
    function polyStr = polyToHTMLString(coeff)
        deg = length(coeff) - 1;  % Determine the degree of the polynomial
        terms = {};                % Cell array to store individual terms

        % Iterate over all coefficients
        for i = 1:length(coeff)
            coef = coeff(i);
            exp = deg - (i - 1);   % Determine the exponent for the current term

            % Skip zero coefficients as they don't affect the expression
            if coef == 0
                continue;
            end

            % Handle signs: add " + " if not the first term;
            % For negative coefficients, add " - " and use the absolute value.
            if coef > 0 && ~isempty(terms)
                term = ' + ';
            elseif coef < 0
                term = ' - ';
                coef = abs(coef);  % Use absolute value for display
            else
                term = '';
            end

            % Show the coefficient unless it's 1 and not the constant term (exp==0)
            if coef ~= 1 || exp == 0
                term = strcat(term, num2str(coef));
            end

            % Add "s" and superscript for exponents greater than 0.
            if exp > 1
                term = strcat(term, 's<sup>', num2str(exp), '</sup>');
            elseif exp == 1
                term = strcat(term, 's');
            end

            % Add the formatted term to the cell array.
            terms{end + 1} = term;
        end

        % Join all terms into a single string
        polyStr = strjoin(terms, '');
        
        % If all coefficients are zero, return "0"
        if isempty(polyStr)
            polyStr = '0';
        end
    end

    % Helper function: Returns an HTML string with the symbolic (variable) preview
  
    function htmlStr = getSymbolicPreviewHTML(ctrlType)
        switch ctrlType
            case 'P'
                numHTML = 'Kp';
                denHTML = '1';
            case 'PI'
                numHTML = 'Kp s + Ki';
                denHTML = 's';
            case 'PD'
                numHTML = 'Kp + Kd s';
                denHTML = '1';
            case 'PID'
                numHTML = 'Kd s<sup>2</sup> + Kp s + Ki';
                denHTML = 's';
            case 'PT1'
                numHTML = 'Kp';
                denHTML = 'T s + 1';
            case 'PIT1'
                numHTML = 'Kp Ki s + Kp';
                denHTML = 'T s<sup>2</sup> + s';
            case 'I2'
                numHTML = 'Ki';
                denHTML = 's<sup>2</sup>';
            case 'PIDT1'
                numHTML = 'Kd s<sup>2</sup> + Kp s + Ki';
                denHTML = 'T s<sup>2</sup> + s';
            case 'Custom'
                numHTML = 'numerator';
                denHTML = 'denominator';
            otherwise
                numHTML = '?';
                denHTML = '?';
        end
        htmlStr = ['<html><head><style>',...
            '.container { display: inline-flex; align-items: center; justify-content: center; width: 100%; }',...
            '.label { font-size:16px; margin-right:10px; }',...
            '.fraction { display:inline-block; text-align:center; }',...
            '.fraction .num { display:block; border-bottom:1px solid black; padding-bottom:2px; }',...
            '.fraction .den { display:block; padding-top:2px; }',...
            '</style></head><body>',...
            '<div class="container">',...
            '<div class="label">Preview:</div>',...
            '<div class="fraction">',...
            '<span class="num">', numHTML, '</span>',...
            '<span class="den">', denHTML, '</span>',...
            '</div>',...
            '</div>',...
            '</body></html>'];
    end
end