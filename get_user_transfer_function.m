function G = get_user_transfer_function()
    % get_user_transfer_function
    % ----------------------------
    % This function creates a modal UI window that allows the user to input
    % the numerator and denominator coefficients for a transfer function.
    % The window contains input fields, a Preview button, a Confirm button,
    % and a Cancel button. When the user clicks Preview, the transfer function
    % is rendered in an HTML preview area on the same line: the text "Preview:"
    % is left-aligned and the polynomial fraction is shown next to it.
    % Decimal commas are allowed (converted to periods). When Confirm is clicked,
    % the transfer function is created and returned.
    
    % Create a modal UI figure
    fig = uifigure('Name', 'Transfer Function Input', ...
                   'Position', [100 100 500 350], 'WindowStyle', 'modal');
               
    % Create a UIHTML component for the preview at the top of the window
    previewHTML = uihtml(fig, 'Position', [20 280 460 50], ...
        'HTMLSource', '<html><body style="font-size:16px; text-align:center;">Preview: </body></html>');
    
    % Create labels and input fields for numerator and denominator
    lblNum = uilabel(fig, 'Position', [20 220 250 22], ...
                     'Text', 'Enter numerator coefficients [b0 b1 ...]:');
    editNum = uieditfield(fig, 'text', 'Position', [20 190 460 22]);
    
    lblDen = uilabel(fig, 'Position', [20 150 250 22], ...
                     'Text', 'Enter denominator coefficients [a0 a1 ...]:');
    editDen = uieditfield(fig, 'text', 'Position', [20 120 460 22]);
    
    % Create three centered buttons at the bottom
    % Total width with buttons and spacing: 3*100 + 2*20 = 340. Left margin = (500-340)/2 = 80.
    btnPreview = uibutton(fig, 'push', 'Text', 'Preview', ...
                          'Position', [80 60 100 30], ...
                          'ButtonPushedFcn', @(btn,event) previewCallback());
    btnConfirm = uibutton(fig, 'push', 'Text', 'Confirm', ...
                          'Position', [200 60 100 30], ...
                          'ButtonPushedFcn', @(btn,event) confirmCallback());
    btnCancel = uibutton(fig, 'push', 'Text', 'Cancel', ...
                         'Position', [320 60 100 30], ...
                         'ButtonPushedFcn', @(btn,event) cancelCallback());
    
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
