function G = get_user_transfer_function()
    % GET_USER_TRANSFER_FUNCTION
    % ----------------------------
    % This function creates a modal UI window that allows the user to input
    % the numerator and denominator coefficients for a transfer function.
    % The window contains input fields, a Preview button, a Confirm button,
    % and a Cancel button with app-matching color scheme.
    
    % Define app color scheme
    appColors = struct(...
        'background', [1 1 1], ...  % White background
        'primary', [0.3 0.5 0.8], ...  % Blue buttons
        'text', [0.2 0.2 0.2], ...  % Dark text
        'accent', [0.95 0.95 0.97], ...  % Light gray-blue panel
        'warn', [0.8 0.3 0.3]);  % Red for cancel
    
    % Create a modal UI figure
    fig = uifigure('Name', 'Transfer Function Input', ...
                   'Position', [100 100 500 350], 'WindowStyle', 'modal', ...
                   'Color', appColors.background);
               
    % Create a UIHTML component for the preview at the top of the window
    previewHTML = uihtml(fig, 'Position', [20 280 460 50], ...
        'HTMLSource', '<html><body style="font-size:16px; text-align:center;">Preview: </body></html>');
    
    % Create labels and input fields for numerator and denominator
    lblNum = uilabel(fig, 'Position', [20 220 250 22], ...
                     'Text', 'Enter numerator coefficients [b0 b1 ...]:', ...
                     'FontColor', appColors.text);
    editNum = uieditfield(fig, 'text', 'Position', [20 190 460 22], ...
                         'BackgroundColor', appColors.accent);
    
    lblDen = uilabel(fig, 'Position', [20 150 250 22], ...
                     'Text', 'Enter denominator coefficients [a0 a1 ...]:', ...
                     'FontColor', appColors.text);
    editDen = uieditfield(fig, 'text', 'Position', [20 120 460 22], ...
                         'BackgroundColor', appColors.accent);
    
    % Create three centered buttons at the bottom
    btnPreview = uibutton(fig, 'push', 'Text', 'Preview', ...
                          'Position', [80 60 100 30], ...
                          'ButtonPushedFcn', @(btn,event) previewCallback(), ...
                          'BackgroundColor', appColors.primary, ...
                          'FontColor', [1 1 1]);
    btnConfirm = uibutton(fig, 'push', 'Text', 'Confirm', ...
                          'Position', [200 60 100 30], ...
                          'ButtonPushedFcn', @(btn,event) confirmCallback(), ...
                          'BackgroundColor', appColors.primary, ...
                          'FontColor', [1 1 1]);
    btnCancel = uibutton(fig, 'push', 'Text', 'Cancel', ...
                         'Position', [320 60 100 30], ...
                         'ButtonPushedFcn', @(btn,event) cancelCallback(), ...
                         'BackgroundColor', appColors.warn, ...
                         'FontColor', [1 1 1]);
    
    % Initialize output
    G = [];
    
    % Callback for the Preview button: update the HTML preview in the same line
    function previewCallback()
        % Retrieve input from the edit fields
        numStr = editNum.Value;
        denStr = editDen.Value;
        
        % Allow decimal commas by replacing commas with periods
        numStr = strrep(numStr, ',', '.');
        denStr = strrep(denStr, ',', '.');
        
        % Convert input strings to numeric vectors
        numVec = str2num(numStr);  %#ok<ST2NM>
        denVec = str2num(denStr);  %#ok<ST2NM>
        
        if isempty(numVec) || isempty(denVec) || any(isnan(numVec)) || any(isnan(denVec))
            previewHTML.HTMLSource = '<html><body style="font-size:16px; text-align:center;">Preview: Invalid input.</body></html>';
        else
            % Convert the coefficient vectors to an HTML-formatted polynomial string
            numHTML = polyToHTMLString(numVec);
            denHTML = polyToHTMLString(denVec);
            % Create an HTML fraction with "Preview:" and the fraction on one horizontal line
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
        end
    end

    % Callback for the Confirm button: validate the input and create the transfer function
    function confirmCallback()
        numStr = strrep(editNum.Value, ',', '.');
        denStr = strrep(editDen.Value, ',', '.');
        numVec = str2num(numStr);  %#ok<ST2NM>
        denVec = str2num(denStr);  %#ok<ST2NM>
        if isempty(numVec) || isempty(denVec) || any(isnan(numVec)) || any(isnan(denVec))
            uialert(fig, 'Invalid input. Please enter valid coefficient vectors.', 'Error');
        else
            G = tf(numVec, denVec);
            uiresume(fig);
            close(fig);
        end
    end

    % Callback for the Cancel button: close the window without creating a transfer function
    function cancelCallback()
        disp('Operation cancelled by user.');
        G = [];
        uiresume(fig);
        close(fig);
    end

    % Block execution until the figure is closed
    uiwait(fig);
end

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