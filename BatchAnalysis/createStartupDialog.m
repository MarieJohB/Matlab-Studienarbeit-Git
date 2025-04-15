function createStartupDialog()
% Create figure for selection matching the provided screenshot
selectionFig = uifigure('Name', 'Control Loop Analysis Tool', 'Position', [500 300 500 400]);
selectionFig.Color = [1 1 1]; % White background

% Define colors for styling
appColors = struct(...
    'headerBg', [0.2 0.4 0.7], ...     % Blue header
    'guiButtonBg', [0.4 0.6 0.9], ...  % Blue button for GUI
    'batchButtonBg', [0.7 0.3 0.7], ... % Purple button for Batch
    'guiPanelBg', [0.9 0.95 1], ...    % Light blue background for GUI section
    'batchPanelBg', [0.95 0.9 1], ...  % Light purple background for Batch section
    'lightText', [1 1 1]);             % White text for dark backgrounds

% Create a header area directly on the figure (no panel)
% Using a colored rectangle
headerArea = uipanel(selectionFig, 'Position', [0 340 500 60], ...
    'BackgroundColor', appColors.headerBg, 'BorderType', 'none');

% Add title to header
uilabel(headerArea, 'Position', [0 0 500 60], 'Text', 'Control Loop Analysis Tool', ...
    'FontSize', 24, 'FontWeight', 'bold', 'FontColor', appColors.lightText, ...
    'HorizontalAlignment', 'center');

% Mode selection area title (directly on the figure)
uilabel(selectionFig, 'Position', [120 300 260 30], 'Text', 'Select Analysis Mode', ...
    'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');

% GUI Mode section (colored rectangle without panel object)
guiSection = uipanel(selectionFig, 'Position', [50 170 400 120], ...
    'BackgroundColor', appColors.guiPanelBg, 'BorderType', 'line', 'Title', '');

% GUI Mode title
uilabel(guiSection, 'Position', [0 90 400 25], 'Text', 'GUI Mode', ...
    'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');

% GUI Mode description
uilabel(guiSection, 'Position', [20 35 360 55], 'Text', ...
    'Interactive single system analysis with real-time feedback and visualization. Best for designing and fine-tuning individual control systems.', ...
    'FontSize', 11, 'WordWrap', 'on');

% Start GUI button
uibutton(guiSection, 'Text', 'Start GUI Mode', 'Position', [130 10 140 30], ...
    'FontSize', 12, 'BackgroundColor', appColors.guiButtonBg, 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) startGUIMode(selectionFig));

% Batch Mode section (colored rectangle without panel object)
batchSection = uipanel(selectionFig, 'Position', [50 40 400 120], ...
    'BackgroundColor', appColors.batchPanelBg, 'BorderType', 'line', 'Title', '');

% Batch Mode title
uilabel(batchSection, 'Position', [0 90 400 25], 'Text', 'Batch Analysis Mode', ...
    'FontWeight', 'bold', 'FontSize', 16, 'HorizontalAlignment', 'center');

% Batch Mode description
uilabel(batchSection, 'Position', [20 35 360 55], 'Text', ...
    'Parameter sweep analysis for optimal controller tuning. Compare multiple configurations and identify stability regions.', ...
    'FontSize', 11, 'WordWrap', 'on');

% Start Batch button
uibutton(batchSection, 'Text', 'Start Batch Mode', 'Position', [130 10 140 30], ...
    'FontSize', 12, 'BackgroundColor', appColors.batchButtonBg, 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) startBatchMode(selectionFig));

% Version info at the bottom
uilabel(selectionFig, 'Position', [0 10 500 20], 'Text', 'v1.1 - Control Engineering Lab', ...
    'FontSize', 9, 'HorizontalAlignment', 'center', 'FontColor', [0.5 0.5 0.5]);
end

function startGUIMode(selectionFig)
% Close the selection dialog
delete(selectionFig);

% Start the main app normally
app = App_automatische_Analyse_eines_Regelkreises_Kopie();
end

function startBatchMode(selectionFig)
% Close the selection dialog
delete(selectionFig);

% Create and show the batch analysis configuration UI
createBatchConfigUI();
end