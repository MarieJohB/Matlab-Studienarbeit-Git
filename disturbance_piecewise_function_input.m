function disturbance_piecewise_function_input(selection_fig, disturbance_label)
% DISTURBANCE_PIECEWISE_FUNCTION_INPUT - Enhanced UI for piecewise disturbance function input
%
% Parameters:
%   selection_fig - The main selection figure handle
%   disturbance_label - Label for the disturbance (e.g. 'd_1')

% Define unified color scheme
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

% Create UI figure with enhanced styling
piece_fig = uifigure('Name', ['Piecewise ' disturbance_label ' Function Input'], 'Position', [300, 150, 800, 650]);
piece_fig.Color = appColors.background;

% Add title panel with enhanced styling
titlePanel = uipanel(piece_fig, 'Position', [10 600 780 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
titleLabel = uilabel(titlePanel, 'Text', ['Piecewise Function Input (' disturbance_label '(t)) - ', num2str(num_sections), ' Sections'], ...
    'Position', [0 0 780 40], 'FontSize', 16, 'FontWeight', 'bold', ...
    'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Add plot area in a panel
previewPanel = uipanel(piece_fig, 'Title', 'Function Preview', ...
    'Position', [10 400 780 190], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

ax = uiaxes(previewPanel, 'Position', [20, 20, 740, 140]);
xlabel(ax, 'Time (s)');
ylabel(ax, disturbance_label);
title(ax, 'Preview');
grid(ax, 'on');

% Create table panel
tablePanel = uipanel(piece_fig, 'Title', 'Section Definitions', ...
    'Position', [10 150 780 240], 'TitlePosition', 'centertop', ...
    'FontWeight', 'bold', 'FontSize', 14, 'BackgroundColor', appColors.panelBg);

% Create table for sections with input fields
columnNames = {'Start Time', 'End Time', 'Time Steps', disturbance_label '(t)'};
columnTypes = {'numeric', 'numeric', 'numeric', 'char'};
columnEditable = [true, true, true, true];
columnWidth = {90, 90, 90, 400};

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

% Create the table with enhanced styling
pieceTable = uitable(tablePanel, 'Position', [20, 20, 740, 190], ...
    'Data', data, 'ColumnName', columnNames, 'ColumnEditable', columnEditable, ...
    'ColumnWidth', columnWidth, 'FontSize', 12, 'RowName', arrayfun(@(x) ['Section ' num2str(x)], 1:num_sections, 'UniformOutput', false));

% Add help button
helpPanel = uipanel(piece_fig, 'Position', [10 90 780 50], 'BorderType', 'none', 'BackgroundColor', appColors.background);

helpBtn = uibutton(helpPanel, 'push', 'Text', 'Function Input Help', ...
    'Position', [340, 10, 120, 30], ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'FontSize', 12, ...
    'ButtonPushedFcn', @(btn, event) showHelpDialog());

% Button panel
buttonPanel = uipanel(piece_fig, 'Position', [10 20 780 60], 'BorderType', 'none', 'BackgroundColor', appColors.background);

% Preview button
previewButton = uibutton(buttonPanel, 'push', 'Position', [240, 15, 100, 30], ...
    'Text', 'Preview', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonPrimary, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_preview_piecewise_function(ax, pieceTable, disturbance_label));

% Confirm button
confirmButton = uibutton(buttonPanel, 'push', 'Position', [350, 15, 100, 30], ...
    'Text', 'Confirm', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonConfirm, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_confirm_piecewise_function(selection_fig, piece_fig, pieceTable));

% Cancel button
cancelButton = uibutton(buttonPanel, 'push', 'Position', [460, 15, 100, 30], ...
    'Text', 'Cancel', ...
    'FontSize', 12, ...
    'BackgroundColor', appColors.buttonCancel, ...
    'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn, event) disturbance_cancel_input(piece_fig));

% Function to show help dialog
function showHelpDialog()
    helpFig = uifigure('Name', 'Function Input Help', 'Position', [350, 200, 500, 500]);
    helpFig.Color = appColors.background;
    
    % Title panel for help
    helpTitlePanel = uipanel(helpFig, 'Position', [10 450 480 40], ...
        'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    helpTitleLabel = uilabel(helpTitlePanel, 'Text', ['Piecewise ' disturbance_label ' Function Reference'], ...
        'Position', [0 0 480 40], 'FontSize', 16, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Help text area
    helpText = uitextarea(helpFig, 'Position', [20 60 460 380], 'Editable', 'off', 'FontSize', 12);
    helpText.Value = {
        ['Piecewise Disturbance Input Guide (' disturbance_label '(t)):'], 
        '----------------------------------------', 
        'For each section, define:', 
        '• Start Time: Beginning of the time segment', 
        '• End Time: End of the time segment', 
        '• Time Steps: Sampling interval within the segment', 
        ['• ' disturbance_label '(t): Mathematical expression defining the disturbance'], 
        '', 
        'Rules for sections:', 
        '• Each section must have End Time > Start Time', 
        '• The Start Time of section n+1 must equal the End Time of section n', 
        '• Time Steps must be positive', 
        '', 
        'Example expressions for disturbances:', 
        '• No disturbance: "0"', 
        '• Constant: "2" or "-0.5"', 
        '• Step functions: "t>=5" (activates at t=5)', 
        '• Pulse: "(t>=2 && t<=3)*5" (pulse of amplitude 5 between t=2 and t=3)', 
        '• Sinusoidal: "0.5*sin(2*t)" (oscillation with amplitude 0.5)', 
        '• Random noise: "0.1*randn(1)"', 
        '• Ramp: "0.2*(t-5)*(t>=5)" (ramp starting at t=5 with slope 0.2)', 
        '', 
        'Example piecewise disturbance:', 
        'Section 1: t=0 to t=5, ' + disturbance_label + '(t) = 0 (no disturbance)', 
        'Section 2: t=5 to t=8, ' + disturbance_label + '(t) = 2 (constant step)', 
        'Section 3: t=8 to t=15, ' + disturbance_label + '(t) = 2-0.5*(t-8) (decreasing ramp)', 
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

% Get section count with enhanced UI
function count = getSectionCount()
    % Create a modal UI figure with styling
    sectionFig = uifigure('Name', 'Number of Sections', 'Position', [400, 300, 350, 250], 'WindowStyle', 'modal');
    sectionFig.Color = appColors.background;
    
    % Title panel
    secTitlePanel = uipanel(sectionFig, 'Position', [10 200 330 40], 'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');
    secTitleLabel = uilabel(secTitlePanel, 'Text', 'Define Number of Sections', ...
        'Position', [0 0 330 40], 'FontSize', 16, 'FontWeight', 'bold', ...
        'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');
    
    % Instruction text
    uilabel(sectionFig, 'Text', ['Enter the number of sections for ' disturbance_label '(t):'], ...
        'Position', [20 160 310 22], 'FontSize', 12, 'HorizontalAlignment', 'center');
    
    % Input panel
    inputPanel = uipanel(sectionFig, 'Position', [75 100 200 50], 'BackgroundColor', appColors.panelBg, 'BorderType', 'line');
    
    % Number input with spinner
    sectionSpinner = uispinner(inputPanel, 'Position', [50 10 100 30], ...
        'Value', 2, 'Limits', [1 5], 'Step', 1, 'FontSize', 14);
    
    % Button panel
    btnPanel = uipanel(sectionFig, 'Position', [10 20 330 60], 'BorderType', 'none', 'BackgroundColor', appColors.background);
    
    % OK button
    okBtn = uibutton(btnPanel, 'push', 'Text', 'OK', ...
        'Position', [125, 15, 100, 30], ...
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