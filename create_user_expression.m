function create_user_expression()
    
    % Check if the figure already exists and close it
    % neccessary for the case that the function gets recalled
    existingFig = findall(0, 'Type', 'figure', 'Name', 'Expression Builder');
    if ~isempty(existingFig)
        close(existingFig);
    end


    % creating a figure for user interaction
    fig = uifigure('Name', 'Expression Builder', 'Position', [100 100 400 500]);

    % Create a scrollable panel
    % needs to be scrollable in case of multiple terms
    scrollPanel = uipanel(fig, 'Position', [10 10 380 480], 'Scrollable', 'on');


    

    % 1. choose how many terms the expression should contain
    uilabel(scrollPanel, 'Position', [20 430 200 22], 'Text', 'Select number of terms:');
    % input for user to enter a value
    inputField = uieditfield(scrollPanel, 'text', 'Position', [150 430 200 22]);

    % Button to get the answer
    submitButton = uibutton(scrollPanel, 'Position', [200 400 100 22], 'Text', 'Submit', 'ButtonPushedFcn', @(btn, event) process_input());
    
    function process_input()

    userInput = inputField.Value;

    number_str = strrep(userInput, ',', '.'); % checking for "," and replacing with "."
    number = str2double(number_str); % converting the input to a double
    

        % testing whether the conversion was succsessfull
        if isnan(number)
         % if not user has to try again
            uiwait(msgbox('Invalid input. Please enter a valid number.', 'Error','error'));
            inputField.Value = '';
            create_user_expression(); % calling the function again to get valid input
        else 
        
            if  mod(number, 1) ~= 0 || number < 1 || number > 5 
                uiwait(msgbox('Invalid input. Please enter a valid integer between 1 and 5.', 'Error','error'));
                inputField.Value = '';
                create_user_expression(); % calling the function again 
            else

            % Delete the old submit button to make space
            delete(submitButton);
           
            % Create input fields and dropdowns to create terms
            yPos = 390; % starting position
            termInputs = struct(); % terms consist of coefficient and time depending term
                for i = 1:number % creating as many terms as user entered priviously

                    % 1) input the coefficient
                    uilabel(scrollPanel, 'Position', [20 yPos 100 22], 'Text', ['Coefficient ', num2str(i), ':']);
                    termInputs(i).coef = uieditfield(scrollPanel, 'numeric', 'Position', [150 yPos 100 22]);

                    yPos = yPos - 30; % place new field underneath
                    
                    % create dropdown so user can select the term
                    uilabel(scrollPanel, 'Position', [20 yPos 100 22], 'Text', ['Term ', num2str(i), ':']);
                    termInputs(i).term = uidropdown(scrollPanel, 'Items', {'t', 't.^2', 't.^3', 't.^4', 'sin(t)', 'cos(t)', 'exp(t)', 'tan(t)', 'atan(t)', 'sqrt(t)', '1/t'}, 'Position', [150 yPos 100 22]);

                    yPos = yPos - 40;
                end

            % Button to create the expression
            % Positioning at the bottom of field
            uibutton(scrollPanel, 'Position', [20 yPos 150 22], 'Text', 'Create Expression', 'ButtonPushedFcn', @(btn, event) create_expression(termInputs)); 
        
            end
        end
        
    end

    % new function to create the time dependent expression
    function create_expression(termInputs)

        % creating an array of strings 
        % array has 1 row
        % array has columns: number of elements in termInputs
        expressionParts = strings(1, length(termInputs)); 

        for i = 1:length(termInputs) % number of terms
            coef = termInputs(i).coef.Value; % extracting the coefficients from the struct
            term = termInputs(i).term.Value; % extracting the terms
            expressionParts(i) = strcat(num2str(coef), '*', term); % the coefficients and terms that belong together are multiplied
        end

        expression = strjoin(expressionParts, ' + '); % all the parts are added to one expression
        expression = char(expression);  % Convert to character array (string scalar)
        target_value = str2func(['@(t)' expression]); % creating a functions that is depending on t
   
        disp(['Created Expression: ', expression]);

        % closing the field at the end
            existingFig = findall(0, 'Type', 'figure', 'Name', 'Expression Builder');
             if ~isempty(existingFig)
             close(existingFig);
             end
        
        % for testing prupose: 
        t = 0:0.01:10;
        y = target_value(t);
        plot(t, y);

    end
end
