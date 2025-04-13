function batchVisualization(batchResults, resultsFile)
% BATCHVISUALIZATION - Creates a visualization UI for batch analysis results
%
% This function creates a tabbed interface showing various visualizations of
% the batch analysis results.
%
% Inputs:
%   batchResults - Structure containing batch analysis results
%   resultsFile - Path to the file where results are saved

% Create figure for results visualization with full HD size
resultsFig = uifigure('Name', 'Batch Analysis Results', 'Position', [0 0 1920 1080]);
resultsFig.Color = [1 1 1];

% Define consistent colors for styling
appColors = struct(...
    'background', [0.95 0.95 0.97], ...     % Light gray background
    'panelHeader', [0.2 0.4 0.7], ...       % Blue header
    'panelBg', [0.95 0.95 0.97], ...        % Light gray panel
    'primary', [0.3 0.6 0.9], ...           % Blue buttons
    'confirm', [0.3 0.8 0.3], ...           % Green confirm button
    'cancel', [0.8 0.3 0.3], ...            % Red cancel button
    'text', [0.2 0.2 0.2], ...              % Dark text
    'lightText', [1 1 1]);                  % White text for dark backgrounds

% Create header panel
headerPanel = uipanel(resultsFig, 'Position', [10 1010 1900 60], ...
    'BackgroundColor', appColors.panelHeader, 'BorderType', 'none');

% Add file info and title to header - centered alignment
[~, fileName, ext] = fileparts(resultsFile);

% Center-aligned title
uilabel(headerPanel, 'Position', [10 30 1880 25], 'Text', 'Batch Analysis Results', ...
    'FontSize', 18, 'FontWeight', 'bold', 'FontColor', appColors.lightText, ...
    'HorizontalAlignment', 'center');

% Center-aligned filename
uilabel(headerPanel, 'Position', [10 5 1880 22], 'Text', ['File: ' fileName ext], ...
    'FontSize', 12, 'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Add Load button directly on the header panel
loadButton = uibutton(headerPanel, 'push', 'Text', 'Load Results', ...
    'Position', [1770 15 120 30], 'FontSize', 12, ...
    'BackgroundColor', [0.1 0.3 0.6], 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) loadResults(resultsFig));

% Create tabs for different analyses with full HD size
resultsTabs = uitabgroup(resultsFig, 'Position', [10 10 1900 990]);

% Summary tab
summaryTab = uitab(resultsTabs, 'Title', 'Summary');
createSummaryTab(summaryTab, batchResults);

% Create tabs for each analysis type
if isfield(batchResults, 'stability')
    stabilityTab = uitab(resultsTabs, 'Title', 'Stability');
    createStabilityTab(stabilityTab, batchResults);
end

if isfield(batchResults, 'nyquist') && ~isempty(batchResults.nyquist{1})
    nyquistTab = uitab(resultsTabs, 'Title', 'Nyquist');
    createNyquistTab(nyquistTab, batchResults);
end

if isfield(batchResults, 'bode') && ~isempty(batchResults.bode{1})
    bodeTab = uitab(resultsTabs, 'Title', 'Bode');
    createBodeTab(bodeTab, batchResults);
end

if isfield(batchResults, 'margins') && ~isempty(batchResults.margins{1})
    marginsTab = uitab(resultsTabs, 'Title', 'Margins');
    createMarginsTab(marginsTab, batchResults);
end

if isfield(batchResults, 'keyParams') && ~isempty(batchResults.keyParams{1})
    keyParamsTab = uitab(resultsTabs, 'Title', 'Key Parameters');
    createKeyParamsTab(keyParamsTab, batchResults);
end

if isfield(batchResults, 'jump') && ~isempty(batchResults.jump{1})
    jumpTab = uitab(resultsTabs, 'Title', 'Jump Analysis');
    createJumpTab(jumpTab, batchResults);
end

% Function to handle loading new results
function loadResults(fig)
    % Prompt for file
    [file, path] = uigetfile('*.mat', 'Load Batch Analysis Results');
    if isequal(file, 0) || isequal(path, 0)
        % User cancelled
        return;
    end
    
    % Full path
    fullPath = fullfile(path, file);
    
    % Load results
    try
        data = load(fullPath);
        if isfield(data, 'batchResults')
            % Ask user what to do
            choice = uiconfirm(fig, 'Do you want to replace the current results or open in a new window?', ...
                'Load Options', 'Options', {'Replace Current', 'New Window', 'Cancel'}, ...
                'DefaultOption', 2, 'CancelOption', 3);
            
            switch choice
                case 'Replace Current'
                    % Close current figure and create new one with the loaded data
                    delete(fig);
                    batchVisualization(data.batchResults, fullPath);
                case 'New Window'
                    % Create a new instance with the loaded data
                    batchVisualization(data.batchResults, fullPath);
                otherwise
                    % Do nothing (cancel)
            end
        else
            uialert(fig, 'The selected file does not contain valid batch analysis results.', 'Invalid File');
        end
    catch ME
        uialert(fig, ['Error loading file: ' ME.message], 'Load Error');
    end
end
end