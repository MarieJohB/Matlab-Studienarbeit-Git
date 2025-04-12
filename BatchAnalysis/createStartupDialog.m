function createStartupDialog()
% Create figure for selection
selectionFig = uifigure('Name', 'Control Loop Analysis Mode', 'Position', [500 500 400 250]);
selectionFig.Color = [1 1 1];

% Find the path to locate images
pathToMLAPP = fileparts(mfilename('fullpath'));

% Try to add the app image (if available)
try
    uiimage(selectionFig, 'Position', [100 180 200 50], ...
        'ImageSource', fullfile(pathToMLAPP, 'App Bilder', 'Standard Regelkreis.png'));
catch
    % If image not found, just continue without it
    disp('App image not found - continuing without it');
end

uilabel(selectionFig, 'Position', [50 150 300 22], 'HorizontalAlignment', 'center', ...
    'Text', 'Select Analysis Mode', 'FontSize', 16, 'FontWeight', 'bold');

% Add buttons for GUI or Batch mode
uibutton(selectionFig, 'Text', 'GUI Mode', 'Position', [100 80 200 40], ...
    'FontSize', 14, 'BackgroundColor', [0.8 0.9 1], ...
    'ButtonPushedFcn', @(btn,event) startGUIMode(selectionFig));

uibutton(selectionFig, 'Text', 'Batch Analysis', 'Position', [100 30 200 40], ...
    'FontSize', 14, 'BackgroundColor', [0.9 0.8 1], ...
    'ButtonPushedFcn', @(btn,event) startBatchMode(selectionFig));
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