function K = get_user_controller()
% get_user_controller
% ---------------------
% This function creates a MATLAB App for selecting a controller type and
% entering its parameters with a live HTML preview. The input dialog layout,
% including the preview area, labels, input fields, and buttons, is modeled
% exactly after the get_user_transfer_function function.
%
% Before any Preview button is pressed, the preview area displays a symbolic
% representation (with variables like Kp, Ki, T, etc.). When Preview is pressed,
% if numeric values are entered, the preview is updated to show the transfer
% function with the actual numbers.
%
% Additionally, if any input field contains a space, the string is truncated
% at the first whitespace so that only the number remains.
%
% The transfer function is created using tf only (no pid objects) and supports
% decimal commas (converted to periods). If the user cancels, the controller
% selection window remains open.
%
% \ac{} \cite{}

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
    
    % Block execution until a valid controller is confirmed
    uiwait(mainFig);
    
    %% Nested Functions
    
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
            catch
                previewHTML.HTMLSource = '<html><body style="font-size:16px; text-align:center;">Preview: Invalid input.</body></html>';
            end
        end
        
        % Callback for Confirm: validate input, create transfer function and close dialog & main figure
        function confirmCallbackK()
            try
                K = computeController(ctrlType, fields);
                uiresume(dlg);
                close(dlg);
                close(mainFig);
            catch
                uialert(dlg, 'Invalid input. Please enter valid parameter values.', 'Error');
            end
        end
        
        % Callback for Cancel: close the dialog and return to the controller selection window
        function cancelCallbackK()
            uiresume(dlg);
            close(dlg);
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
        deg = length(coeff) - 1;
        polyStr = '';
        for i = 1:length(coeff)
            coef = coeff(i);
            if coef == 0
                continue;
            end
            if ~isempty(polyStr)
                if coef > 0
                    polyStr = [polyStr, ' + '];  %#ok<AGROW>
                else
                    polyStr = [polyStr, ' - '];  %#ok<AGROW>
                    coef = abs(coef);
                end
            end
            expVal = deg - i + 1;
            if expVal > 1
                polyStr = [polyStr, num2str(coef), ' s<sup>', num2str(expVal), '</sup>'];
            elseif expVal == 1
                polyStr = [polyStr, num2str(coef), ' s'];
            else
                polyStr = [polyStr, num2str(coef)];
            end
        end
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
