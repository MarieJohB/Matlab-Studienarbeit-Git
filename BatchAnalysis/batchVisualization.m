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
fileLabel = uilabel(headerPanel, 'Position', [10 5 1880 22], 'Text', ['File: ' fileName ext], ...
    'FontSize', 12, 'FontColor', appColors.lightText, 'HorizontalAlignment', 'center');

% Adjust Load button position to make room for New Analysis button
loadButton = uibutton(headerPanel, 'push', 'Text', 'Load Results', ...
    'Position', [1640 15 120 30], 'FontSize', 12, ...
    'BackgroundColor', [0.1 0.3 0.6], 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) loadResults());

% Add New Analysis button next to Load button
newAnalysisButton = uibutton(headerPanel, 'push', 'Text', 'New Analysis', ...
    'Position', [1770 15 120 30], 'FontSize', 12, ...
    'BackgroundColor', [0.1 0.3 0.6], 'FontColor', appColors.lightText, ...
    'ButtonPushedFcn', @(btn,event) startNewAnalysis());

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
function loadResults()
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
            % Store a reference to the current figure
            currentFig = resultsFig;
            
            % Ask user what to do
            choice = uiconfirm(currentFig, 'Do you want to replace the current results or open in a new window?', ...
                'Load Options', 'Options', {'Replace Current', 'New Window', 'Cancel'}, ...
                'DefaultOption', 2, 'CancelOption', 3);
            
            switch choice
                case 'Replace Current'
                    % Update the current figure with the loaded data
                    % First clear all tabs
                    delete(resultsTabs.Children);
                    
                    % Update file info in header
                    [~, newFileName, newExt] = fileparts(fullPath);
                    fileLabel.Text = ['File: ' newFileName newExt];
                    
                    % Recreate tabs for the new data
                    resultsTabs = uitabgroup(resultsFig, 'Position', [10 10 1900 990]);
                    
                    % Summary tab
                    summaryTab = uitab(resultsTabs, 'Title', 'Summary');
                    createSummaryTab(summaryTab, data.batchResults);
                    
                    % Create tabs for each analysis type
                    if isfield(data.batchResults, 'stability')
                        stabilityTab = uitab(resultsTabs, 'Title', 'Stability');
                        createStabilityTab(stabilityTab, data.batchResults);
                    end
                    
                    if isfield(data.batchResults, 'nyquist') && ~isempty(data.batchResults.nyquist{1})
                        nyquistTab = uitab(resultsTabs, 'Title', 'Nyquist');
                        createNyquistTab(nyquistTab, data.batchResults);
                    end
                    
                    if isfield(data.batchResults, 'bode') && ~isempty(data.batchResults.bode{1})
                        bodeTab = uitab(resultsTabs, 'Title', 'Bode');
                        createBodeTab(bodeTab, data.batchResults);
                    end
                    
                    if isfield(data.batchResults, 'margins') && ~isempty(data.batchResults.margins{1})
                        marginsTab = uitab(resultsTabs, 'Title', 'Margins');
                        createMarginsTab(marginsTab, data.batchResults);
                    end
                    
                    if isfield(data.batchResults, 'keyParams') && ~isempty(data.batchResults.keyParams{1})
                        keyParamsTab = uitab(resultsTabs, 'Title', 'Key Parameters');
                        createKeyParamsTab(keyParamsTab, data.batchResults);
                    end
                    
                    if isfield(data.batchResults, 'jump') && ~isempty(data.batchResults.jump{1})
                        jumpTab = uitab(resultsTabs, 'Title', 'Jump Analysis');
                        createJumpTab(jumpTab, data.batchResults);
                    end
                    
                case 'New Window'
                    % Open in a new window without closing the current one
                    % Create a new instance with separate scoping to prevent closure of current window
                    openInNewWindow(data.batchResults, fullPath);
                    
                otherwise
                    % Do nothing (cancel)
            end
        else
            uialert(resultsFig, 'The selected file does not contain valid batch analysis results.', 'Invalid File');
        end
    catch ME
        uialert(resultsFig, ['Error loading file: ' ME.message], 'Load Error');
    end
end

% Function to start a new batch analysis - does NOT close current window
function startNewAnalysis()
    % Start new batch configuration without closing current window
    createBatchConfigUI();
end

end

% External helper function to open a new window without closing the current one
function openInNewWindow(batchResults, resultsFile)
    % This function creates a new, separate instance of the visualization
    % It is kept as a separate function to ensure the calling context is isolated
    
    % Call batchVisualization in a completely separate context
    batchVisualization(batchResults, resultsFile);
end