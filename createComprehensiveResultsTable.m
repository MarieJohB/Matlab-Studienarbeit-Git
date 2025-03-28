function htmlContent = createComprehensiveResultsTable(successCount, numSimulations, metricsData, correlationData)
    % CREATECOMPREHENSIVERESULTSTABLE - Create comprehensive HTML table for all metrics
    % 
    % Enhanced version with horizontal layout for statistics to improve readability
    % while maintaining the same table width as other tables
    %
    % Parameters:
    %   successCount - Number of successful simulations
    %   numSimulations - Total number of simulations
    %   metricsData - Cell array of metric data with format:
    %                 {MetricName, Mean, StdDev, CoV, Unit, Min, Max}
    %   correlationData - Cell array of correlation data with format:
    %                     {MetricName, KpCorr, KiCorr, KdCorr}
    %
    % Returns:
    %   htmlContent - HTML string for the results visualization
    
    % HTML Header with styling - matching other tables in the app
    htmlHeader = [...
        '<html><head><style>', ...
        'body {font-family: Arial, sans-serif; font-size: 12px; margin: 0; padding: 10px; color: #333;}', ...
        'h3 {color: #446699; margin-top: 0; margin-bottom: 10px;}', ...
        'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; margin-bottom: 15px; }', ...
        'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 14px; }', ...
        'td { border: 1px solid #DDDDDD; padding: 8px; font-size: 13px; }', ...
        'td:first-child { font-weight: bold; }', ...
        'tr:nth-child(even) { background-color: #F2F2F2; }', ...
        'tr:hover { background-color: #E8F4F8; }', ...
        '.section-header { background-color: #5D87B1; color: white; font-weight: bold; }', ...
        '.success-rate { color: #2ECC71; font-weight: bold; }', ...
        '.failure-rate { color: #E74C3C; font-weight: bold; }', ...
        '.correlation-positive { color: #27AE60; }', ...
        '.correlation-negative { color: #C0392B; }', ...
        '.high-variability { color: #E74C3C; font-weight: bold; }', ...
        '.low-variability { color: #2ECC71; font-weight: bold; }', ...
        '.metric-header { background-color: #4472C4; color: white; font-weight: bold; }', ...
        '</style></head><body>'];
        
    % Main heading
    htmlBody = '<h3>Monte Carlo Simulation Results: All Metrics</h3>';
    
    % Simulation info with success rate
    successRate = 100 * successCount / numSimulations;
    htmlBody = [htmlBody, ...
        '<table>', ...
        '<tr><th colspan="2">Simulation Information</th></tr>', ...
        '<tr><td width="50%">Total Simulations</td><td>', num2str(numSimulations), '</td></tr>', ...
        '<tr><td>Successful Simulations</td><td>', num2str(successCount), ...
        ' <span class="success-rate">(', num2str(successRate, '%.1f'), '%)</span></td></tr>'];
    
    if successCount < numSimulations
        failureRate = 100 * (numSimulations - successCount) / numSimulations;
        htmlBody = [htmlBody, ...
            '<tr><td>Failed Simulations</td><td>', num2str(numSimulations - successCount), ...
            ' <span class="failure-rate">(', num2str(failureRate, '%.1f'), '%)</span></td></tr>'];
    end
    
    htmlBody = [htmlBody, '</table>'];
    
    % Create horizontal metrics table - new layout with metrics as columns
    htmlBody = [htmlBody, ...
        '<table>', ...
        '<tr><th colspan="', num2str(length(metricsData) + 1), '">Performance Metrics Summary</th></tr>', ...
        '<tr class="metric-header">', ...
        '<td style="color: black;">Statistic</td>'];
    
    % Add metric names as column headers with explicit black color
    for i = 1:length(metricsData)
        htmlBody = [htmlBody, '<td style="color: black;">', metricsData{i}{1}, '</td>'];
    end
    
    htmlBody = [htmlBody, '</tr>'];
    
    % Define the statistics to show
    statTypes = {'Mean', 'Std Dev', 'CoV', 'Min', 'Max'};
    
    % For each statistic type, create a row with values for all metrics
    for statIdx = 1:length(statTypes)
        stat = statTypes{statIdx};
        htmlBody = [htmlBody, '<tr><td style="color: black;">', stat, '</td>'];
        
        % Add values for each metric
        for i = 1:length(metricsData)
            metric = metricsData{i};
            unit = metric{5}; % Unit is stored in position 5
            
            % Based on stat type, pick the right column from metrics data
            switch stat
                case 'Mean'
                    value = metric{2}; % Mean is in position 2
                    htmlBody = [htmlBody, '<td>', value, ' ', unit, '</td>'];
                case 'Std Dev'
                    value = metric{3}; % Std Dev is in position 3
                    htmlBody = [htmlBody, '<td>', value, ' ', unit, '</td>'];
                case 'CoV'
                    value = metric{4}; % CoV is in position 4
                    % Add special styling for CoV
                    covValue = str2double(strrep(value, '%', ''));
                    if isnan(covValue)
                        covClass = '';
                    elseif covValue > 25
                        covClass = 'high-variability';
                    else
                        covClass = 'low-variability';
                    end
                    htmlBody = [htmlBody, '<td class="', covClass, '">', value, '</td>'];
                case 'Min'
                    value = metric{6}; % Min is in position 6
                    htmlBody = [htmlBody, '<td>', value, ' ', unit, '</td>'];
                case 'Max'
                    value = metric{7}; % Max is in position 7
                    htmlBody = [htmlBody, '<td>', value, ' ', unit, '</td>'];
            end
        end
        
        htmlBody = [htmlBody, '</tr>'];
    end
    
    htmlBody = [htmlBody, '</table>'];
    
    % Create horizontal correlation table - improved layout with parameters as rows
    htmlBody = [htmlBody, ...
        '<table>', ...
        '<tr><th colspan="', num2str(length(correlationData) + 1), '">Parameter Correlations</th></tr>', ...
        '<tr class="metric-header">', ...
        '<td style="color: black;">Parameter</td>'];
    
    % Add metric names as column headers with explicit black color
    for i = 1:length(correlationData)
        htmlBody = [htmlBody, '<td style="color: black;">', correlationData{i}{1}, '</td>'];
    end
    
    htmlBody = [htmlBody, '</tr>'];
    
    % Define the parameters
    params = {'Kp', 'Ki', 'Kd'};
    
    % For each parameter, create a row with correlation values for all metrics
    for paramIdx = 1:length(params)
        param = params{paramIdx};
        htmlBody = [htmlBody, '<tr><td style="color: black;">', param, '</td>'];
        
        % Add correlation values for each metric
        for i = 1:length(correlationData)
            corr = correlationData{i};
            
            % Get correlation value (paramIdx + 1 because first element is metric name)
            value = corr{paramIdx + 1};
            
            % Determine color class based on correlation value
            corrValue = str2double(value);
            if isnan(corrValue)
                colorClass = '';
            elseif corrValue > 0
                colorClass = 'correlation-positive';
            else
                colorClass = 'correlation-negative';
            end
            
            % Add cell with correlation value
            htmlBody = [htmlBody, '<td class="', colorClass, '">', value, '</td>'];
        end
        
        htmlBody = [htmlBody, '</tr>'];
    end
    
    htmlBody = [htmlBody, '</table>'];
    
    % Interpretation notes - keep as-is but ensure width matches
    htmlBody = [htmlBody, ...
        '<table>', ...
        '<tr><th>Interpretation Notes</th></tr>', ...
        '<tr><td>', ...
        'Correlation interpretation:<br>', ...
        '• <span class="correlation-positive">Positive values</span>: Parameter increases → metric increases<br>', ...
        '• <span class="correlation-negative">Negative values</span>: Parameter decreases → metric decreases<br>', ...
        '• Values close to ±1 indicate stronger relationships<br>', ...
        '• Values close to 0 indicate weak or no relationship<br><br>', ...
        'Coefficient of Variation (CoV):<br>', ...
        '• <span class="low-variability">Low values (<25%)</span>: More consistent/predictable results<br>', ...
        '• <span class="high-variability">High values (>25%)</span>: More variable/unpredictable results', ...
        '</td></tr>', ...
        '</table>'];
    
    % Complete HTML
    htmlContent = [htmlHeader, htmlBody, '</body></html>'];