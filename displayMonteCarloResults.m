function htmlContent = displayMonteCarloResults(metric, successCount, numSimulations, metricsData, correlationsData, metricNames, showAllMetrics)
    % DISPLAYMONTECARLOESULTS - Display Monte Carlo simulation results
    %
    % This function handles the display of Monte Carlo simulation results,
    % supporting both single-metric and all-metrics views.
    %
    % Inputs:
    %   metric - Current metric name for single-metric view
    %   successCount - Number of successful simulations
    %   numSimulations - Total number of simulations
    %   metricsData - Cell array of metric statistics 
    %   correlationsData - Cell array of correlation data
    %   metricNames - Cell array of all available metric names (optional)
    %   showAllMetrics - Boolean flag to show all metrics (default: false)
    %
    % Output:
    %   htmlContent - HTML string containing the formatted results
    
    % Default to single metric view if not specified
    if nargin < 7
        showAllMetrics = false;
    end
    
    % Default metric names if not provided
    if nargin < 6 || isempty(metricNames)
        metricNames = {metric};
    end
    
    try
        if showAllMetrics
            % Transform data to comprehensive format if needed
            [transformedMetricsData, transformedCorrelationsData] = ...
                transformMonteCarloData(metricsData, correlationsData, metricNames);
            
            % Use comprehensive results table
            htmlContent = createComprehensiveResultsTable(successCount, numSimulations, ...
                transformedMetricsData, transformedCorrelationsData);
        else
            % Single metric view - use standard results table
            htmlContent = createResultsTable(metric, successCount, numSimulations, ...
                metricsData, correlationsData);
        end
    catch ME
        % Handle errors with informative message
        errorMsg = ['<html><body style="font-family: Arial, sans-serif;">', ...
            '<h3 style="color: #d9534f;">Error Displaying Results</h3>', ...
            '<p>An error occurred while displaying the Monte Carlo simulation results:</p>', ...
            '<pre style="background-color: #f8f8f8; padding: 10px; border: 1px solid #ddd;">', ...
            ME.message, '<br>', ...
            ME.stack(1).name, ' (line ', num2str(ME.stack(1).line), ')', ...
            '</pre>', ...
            '<p>Please check the console for more details.</p>', ...
            '</body></html>'];
        
        % Display detailed error in console
        disp('Error in displayMonteCarloResults:');
        disp(getReport(ME, 'extended'));
        
        htmlContent = errorMsg;
    end
end