function batchResults = runBatchAnalysisFromCommandLine(G, K, paramInfo, analysisOptions, savePath)
    % Create a headless app instance
    app = App_automatische_Analyse_eines_Regelkreises_Kopie();
    app.UIFigure.Visible = 'off';
    
    % Run batch analysis without progress UI
    disp('Running batch analysis...');
    disp(['Parameter: ' paramInfo.type ' ' paramInfo.coeffType ' [' num2str(paramInfo.index) ']']);
    disp(['Range: ' num2str(paramInfo.min) ' to ' num2str(paramInfo.max) ' with step ' num2str(paramInfo.step)]);
    
    % Run analysis
    batchResults = runBatchAnalysis(app, G, K, paramInfo, analysisOptions, savePath);
    
    % Clean up app
    delete(app);
    
    disp('Analysis complete!');
    disp(['Results saved to: ' savePath]);
    
    % Return results
    return;
end