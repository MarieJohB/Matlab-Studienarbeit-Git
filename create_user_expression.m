function create_user_expression()
    
    % Check if the figure already exists and close it
    % neccessary for the case that the function gets recalled
    existingFig = findall(0, 'Type', 'figure', 'Name', 'Expression Builder');
    if ~isempty(existingFig)
        close(existingFig);
    end


    % creating a figure for user interaction
    fig = uifigure('Name', 'Expression Builder');
    
    % 1. choose how many terms the expression should contain
    uilabel(fig, 'Position', [20 350 200 22], 'Text', 'Select number of terms:');
    % input for user to enter a value
    inputField = uieditfield(fig, 'text', 'Position', [150 350 200 22]);

    % Button to get the answer
    submitButton = uibutton(fig, 'Position', [200 300 100 22], 'Text', 'Submit', 'ButtonPushedFcn', @(btn, event) process_input());
    
    function process_input()

    userInput = inputField.Value;

    number_str = strrep(userInput, ',', '.'); % checking for "," and replacing with "."
    number = str2double(number_str); % converting the input to a double
    

        % testing whether the conversion was succsessfull
        if isnan(number)
         % if not user has to try again
            uiwait(msgbox('Invalid input. Please enter a valid number.', 'Error','error'));
            inputField.Value = '';
            create_user_expression(); % calling the function again 
        else 
        
            if  mod(number, 1) ~= 0 || number < 1 || number > 10 
                uiwait(msgbox('Invalid input. Please enter a valid integer between 1 and 10.', 'Error','error'));
                inputField.Value = '';
                create_user_expression(); % calling the function again 
            else

            % Delete the old submit button to make space
            delete(submitButton);
           
            % Create input fields and dropdowns to create terms
            yPos = 300;
            termInputs = struct();
                for i = 1:number
                    uilabel(fig, 'Position', [20 yPos 100 22], 'Text', ['Coefficient ', num2str(i), ':']);
                    termInputs(i).coef = uieditfield(fig, 'numeric', 'Position', [150 yPos 100 22]);

                    yPos = yPos - 40;

                    uilabel(fig, 'Position', [20 yPos 100 22], 'Text', ['Term ', num2str(i), ':']);
                    termInputs(i).term = uidropdown(fig, 'Items', {'t', 't^2', 'sin(t)', 'cos(t)', 'exp(t)'}, 'Position', [150 yPos 100 22]);

                    yPos = yPos - 40;
                end

            % Button to create the expression
            uibutton(fig, 'Position', [390 50 120 22], 'Text', 'Create Expression', 'ButtonPushedFcn', @(btn, event) create_expression(termInputs));
            end
        end
    end

    function create_expression(termInputs)


        expressionParts = strings(1, length(termInputs));

        for i = 1:length(termInputs)
            coef = termInputs(i).coef.Value;
            term = termInputs(i).term.Value;
            expressionParts(i) = strcat(num2str(coef), '*', term);
        end

        expression = strjoin(expressionParts, ' + ');
        expression = char(expression);  % Convert to character array (string scalar)
        target_value = str2func(['@(t)' expression]);
        disp(['Created Expression: ', expression]);
        % Optional: Further processing with target_value
    end
end
