function batchVisualization(batchResults, resultsFile)
% BATCHVISUALIZATION - Creates a visualization UI for batch analysis results
%
% This function creates a tabbed interface showing various visualizations of
% the batch analysis results.
%
% Inputs:
%   batchResults - Structure containing batch analysis results
%   resultsFile - Path to the file where results are saved

% Create figure for results visualization
resultsFig = uifigure('Name', 'Batch Analysis Results', 'Position', [100 100 1000 700]);
resultsFig.Color = [1 1 1];

% Add file info
[~, fileName, ext] = fileparts(resultsFile);
uilabel(resultsFig, 'Position', [20 670 960 20], 'Text', ['Results File: ' fileName ext], ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

% Create tabs for different analyses
resultsTabs = uitabgroup(resultsFig, 'Position', [10 10 980 650]);

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
end