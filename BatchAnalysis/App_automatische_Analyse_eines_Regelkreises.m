function App_automatische_Analyse_eines_Regelkreises(startMode)
% This is the main entry point for the app
% It can be called in multiple ways:
%   1. App_automatische_Analyse_eines_Regelkreises() - shows selection dialog
%   2. App_automatische_Analyse_eines_Regelkreises('gui') - starts GUI directly
%   3. App_automatische_Analyse_eines_Regelkreises('batch') - starts batch mode directly

if nargin == 0
    % No arguments - show selection dialog
    createStartupDialog();
elseif nargin == 1 && strcmp(startMode, 'gui')
    % Start in GUI mode
    app = App_automatische_Analyse_eines_Regelkreises_Kopie();
elseif nargin == 1 && strcmp(startMode, 'batch')
    % Start in batch mode
    createBatchConfigUI();
else
    % Invalid arguments
    error('Invalid argument. Use: App_automatische_Analyse_eines_Regelkreises(), App_automatische_Analyse_eines_Regelkreises(''gui''), or App_automatische_Analyse_eines_Regelkreises(''batch'')');
end
end