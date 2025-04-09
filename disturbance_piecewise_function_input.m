function disturbance_piecewise_function_input(selection_fig, disturbance_label)
% DISTURBANCE_PIECEWISE_FUNCTION_INPUT - Enhanced UI for piecewise disturbance function input
% Designed to match the dimensions and layout of target_piecewise_function_input
%
% Parameters:
%   selection_fig - The main selection figure handle
%   disturbance_label - Label for the disturbance (e.g. 'd_1' or 'd_2')

% Define unified color scheme to match continuous function input
appColors = struct(...
    'background', [0.95 0.95 0.97], ...    % Light gray background
    'panelHeader', [0.2 0.4 0.7], ...      % Blue panel header
    'panelBg', [0.95 0.95 0.97], ...       % Light panel background
    'buttonPrimary', [0.3 0.6 0.9], ...    % Blue buttons
    'buttonConfirm', [0.3 0.8 0.3], ...    % Green confirm button
    'buttonCancel', [0.8 0.3 0.3], ...     % Red cancel button
    'buttonHelp', [0.5 0.5 0.5], ...       % Gray help button
    'text', [0.2 0.2 0.2], ...             % Dark text
    'lightText', [1 1 1]);                 % White text for dark backgrounds

% First, get the number of sections with improved UI
num_sections = getSectionCount();

% If user cancelled section selection, exit
if isempty(num_sections)
    disp('Operation cancelled by user.');
    return;
end

% Create UI figure with enhanced styling - SAME SIZE AS TARGET FUNCTION INPUT
piece_fig = uifigure('Name', ['Piecewise ', disturbance_label, ' Function Input'], 'Position', [400, 200, 700, 615]);
piece_fig.Color = appColors.background;

% Add title panel with enhanced styling - MATCH TARGET FUNCTION INPUT
titlePanel = uipanel(piece_fig, 'Position', [10 565 680 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Text', ['Piecewise Function Input (', disturbance_label, '(t)) - ', num2str(num_sections), ' Sections'], ...
    'Position', [0 0 680 40], 'FontSize', 16, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Add plot area in a panel - SMALLER TO MAKE ROOM FOR SECTION DEFINITIONS
previewPanel = uipanel(piece_fig, 'Title', 'Function Preview', ...
    'Position', [10 300 680 255], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

ax = uiaxes(previewPanel, 'Position', [20, 10, 640, 205]);
xlabel(ax, 'Time (s)');
% No y-label as requested
grid(ax, 'on');

% Create table panel - ENLARGED TO USE NEWLY AVAILABLE SPACE
tablePanel = uipanel(piece_fig, 'Title', 'Section Definitions', ...
    'Position', [10 105 680 185], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Create table for sections with input fields
columnNames = {'Start Time', 'End Time', 'Time Steps', [disturbance_label, '(t)']};
columnTypes = {'numeric', 'numeric', 'numeric', 'char'};
columnEditable = [true, true, true, true];
columnWidth = {90, 90, 90, 310};

% Initialize data with default values
data = cell(num_sections, 4);
for i = 1:num_sections
    if i == 1
        data{i, 1} = 0; % Default start time for first section
    else
        data{i, 1} = data{i-1, 2}; % Start time = end time of previous section
    end
    data{i, 2} = data{i, 1} + 10; % Default end time
    data{i, 3} = 0.1; % Default time steps
    data{i, 4} = ''; % Empty function by default
end

% Create the table with enhanced styling - MADE TALLER TO USE AVAILABLE SPACE
pieceTable = uitable(tablePanel, 'Position', [20, 10, 640, 145], ...
    'Data', data, 'ColumnName', columnNames, 'ColumnEditable', columnEditable, ...
    'ColumnWidth', columnWidth, 'FontSize', 12, 'RowName', arrayfun(@(x) ['Section ' num2str(x)], 1:num_sections, 'UniformOutput', false));

% Add cell edit callback to synchronize time steps and update start times
pieceTable.CellEditCallback = @(src, event) updateTableValues(src, event);

% Button panel - MATCH TARGET FUNCTION INPUT
buttonPanel = uipanel(piece_fig, 'Title', 'Actions', ...
    'Position', [10 10 680 85], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Center buttons in panel - MATCH TARGET FUNCTION INPUT
panelCenter = 680/2;
buttonWidth = 100;
spacing = 20;
totalWidth = 3*buttonWidth + 2*spacing;
startX = panelCenter - totalWidth/2;

% Preview button - MATCH TARGET FUNCTION INPUT
previewButton = uibutton(buttonPanel, 'push', 'Position', [startX, 25, buttonWidth, 30], ...
    'Text', 'Preview', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_preview_piecewise_function(ax, pieceTable, disturbance_label));

% Confirm button - MATCH TARGET FUNCTION INPUT
confirmButton = uibutton(buttonPanel, 'push', 'Position', [startX + buttonWidth + spacing, 25, buttonWidth, 30], ...
    'Text', 'Confirm', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonConfirm, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_confirm_piecewise_function(selection_fig, piece_fig, pieceTable));

% Cancel button - MATCH TARGET FUNCTION INPUT
cancelButton = uibutton(buttonPanel, 'push', 'Position', [startX + 2*(buttonWidth + spacing), 25, buttonWidth, 30], ...
    'Text', 'Cancel', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonCancel, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_cancel_input(piece_fig));

% Help button centered under the Confirm button - MATCH TARGET FUNCTION INPUT
helpBtn = uibutton(buttonPanel, 'push', 'Text', 'Help', ...
    'Position', [startX + buttonWidth + spacing, 5, buttonWidth, 18], ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'FontSize', 10, ...
    'ButtonPushedFcn', @(btn, event) showHelpDialog());

% Function to update time steps and start times when cells are edited
function updateTableValues(src, event)
    % Get current table data
    data = src.Data;
    
    % Check which column was edited
    if event.Indices(2) == 3 % Time Steps column
        % Get the new time step value
        newTimeStepStr = event.NewData;
        if ischar(newTimeStepStr)
            newTimeStepStr = strrep(newTimeStepStr, ',', '.');
            newTimeStep = str2double(newTimeStepStr);
        else
            newTimeStep = event.NewData;
        end
        
        % Update all rows with the same time step
        for row = 1:size(data, 1)
            data{row, 3} = newTimeStep;
        end
        
        % Update the table data
        src.Data = data;
    elseif event.Indices(2) == 2 % End Time column
        % If end time is edited, update the start time of the next section
        row = event.Indices(1);
        
        % Convert comma to period if needed
        if ischar(event.NewData)
            newData = str2double(strrep(event.NewData, ',', '.'));
        else
            newData = event.NewData;
        end
        
        if row < size(data, 1)
            data{row+1, 1} = newData; % Set start time of next section
            
            % Also update any subsequent sections to maintain continuity
            for nextRow = row+2:size(data, 1)
                data{nextRow, 1} = data{nextRow-1, 2};
            end
            
            src.Data = data;
        end
    end
end

% Function to show help dialog - Customized for disturbances
    function showHelpDialog()
    helpFig = uifigure('Name', 'Function Input Help', 'Position', [350, 200, 500, 500]);
    helpFig.Color = appColors.background;
    
    % Title panel for help
    helpTitlePanel = uipanel(helpFig, 'Position', [10 450 480 40], ...
        'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    helpTitleLabel = uilabel(helpTitlePanel, 'Text', 'Piecewise Function Reference', ...
        'Position', [0 0 480 40], 'FontSize', 16, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Help text area
    helpText = uitextarea(helpFig, 'Position', [20 60 460 380], 'Editable', 'off', 'FontSize', 12);
    helpText.Value = {
        'Piecewise Function Input Guide:', 
        '----------------------------------------', 
        'For each section, define:', 
        '• Start Time: Beginning of the time segment', 
        '• End Time: End of the time segment', 
        '• Time Steps: Sampling interval within the segment', 
        '• Function: Mathematical expression defining the function', 
        '', 
        'Rules for sections:', 
        '• Each section must have End Time > Start Time', 
        '• The Start Time of section n+1 must equal the End Time of section n', 
        '• Time Steps must be positive', 
        '', 
        'Example expressions:', 
        '• Constant: "5" or "3.14"', 
        '• Linear functions: "2*t" or "t+5"', 
        '• Polynomials: "t^2 + 3*t + 1"', 
        '• Trigonometric: "sin(t)" or "cos(2*t)"', 
        '• Exponentials: "exp(-t)" or "exp(-0.5*t)*sin(t)"', 
        '• Step functions: "t>=2" (returns 1 when t≥2, 0 otherwise)', 
        '• Pulse: "(t>=2 && t<=3)*5" (pulse of amplitude 5 between t=2 and t=3)', 
        '', 
        'Example piecewise function:', 
        'Section 1: t=0 to t=5, function = t^2', 
        'Section 2: t=5 to t=10, function = 25 (constant)', 
        'Section 3: t=10 to t=15, function = 25-5*(t-10) (decreasing line)', 
        '', 
        'Tips:', 
        '• Use "*" for multiplication: "5*t" not "5t"', 
        '• Ensure proper syntax for functions and operators'
    };
    
    % Close button
    closeBtn = uibutton(helpFig, 'push', 'Text', 'Close', ...
        'Position', [200, 20, 100, 30], ...
        'BackgroundColor', appColors.buttonPrimary, ...
        'FontColor', appColors.lightText, ...
        'ButtonPushedFcn', @(~,~) close(helpFig));
end

% Get section count with enhanced UI - Matching target_piecewise_function_input
function count = getSectionCount()
    % Create a modal UI figure with styling matching the signal type selection window
    sectionFig = uifigure('Name', 'Number of Sections', 'Position', [500, 300, 300, 220], 'WindowStyle', 'modal');
    sectionFig.Color = appColors.background;
    
    % Add title panel with enhanced styling - matches signal type selection
    titlePanel = uipanel(sectionFig, 'Position', [10 170 280 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    titleLabel = uilabel(titlePanel, 'Text', 'Select Number of Sections', ...
        'Position', [0 0 280 40], 'FontSize', 16, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Create a larger spinner instead of panel with spinner
    % Position it in the center of the window
    sectionSpinner = uispinner(sectionFig, 'Position', [100, 100, 100, 40], ...
        'Value', 2, 'Limits', [1 5], 'Step', 1, 'FontSize', 14);
    
    % OK button - matching the Cancel button position in signal type window
    okBtn = uibutton(sectionFig, 'push', 'Text', 'OK', ...
        'Position', [100, 10, 100, 30], ...
        'BackgroundColor', appColors.buttonConfirm, ...
        'FontColor', appColors.lightText, ...
        'FontSize', 12, ...
        'ButtonPushedFcn', @(~,~) confirmSections());
    
    % Initialize output
    count = [];
    
    % Handle figure close request
    sectionFig.CloseRequestFcn = @(~,~) closeFigure();
    
    % Confirm button callback
    function confirmSections()
        count = sectionSpinner.Value;
        uiresume(sectionFig);
        delete(sectionFig);
    end
    
    % Close figure callback
    function closeFigure()
        count = [];
        uiresume(sectionFig);
        delete(sectionFig);
    end
    
    % Wait for user response
    uiwait(sectionFig);
end

% Handle window close event
piece_fig.CloseRequestFcn = @(src, event) disturbance_cancel_input(src);

% Wait for user to confirm or cancel
uiwait(piece_fig);
end