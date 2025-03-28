function [transformedMetricsData, transformedCorrelationsData] = transformMonteCarloData(metricData, correlationData, metricNames)
    % TRANSFORMMONTECARLODATA - Transform Monte Carlo simulation data between formats
    %
    % This function transforms metric and correlation data between single-metric
    % and all-metrics formats to ensure compatibility with different display functions.
    %
    % Inputs:
    %   metricData - Original metrics data (can be in various formats)
    %   correlationData - Original correlation data (can be in various formats)
    %   metricNames - Cell array of metric names (optional)
    %
    % Outputs:
    %   transformedMetricsData - Standardized metrics data
    %   transformedCorrelationsData - Standardized correlation data
    
    % Initialize outputs
    transformedMetricsData = {};
    transformedCorrelationsData = {};
    
    % If inputs are empty, return empty outputs
    if isempty(metricData) || ~iscell(metricData)
        return;
    end
    
    % Determine data format by examining the structure
    if ~isempty(metricData) && iscell(metricData{1})
        % Get dimensions of the first element
        firstMetric = metricData{1};
        
        % Check if we have a structure with "all metrics" format
        if length(firstMetric) >= 6 && ischar(firstMetric{1})
            % Already in comprehensive format [name, mean, stddev, unit, min, max]
            transformedMetricsData = metricData;
            
            % Handle correlation data if present
            if ~isempty(correlationData) && iscell(correlationData{1}) && length(correlationData{1}) >= 4
                transformedCorrelationsData = correlationData;
            else
                % Transform correlation data to comprehensive format if needed
                transformedCorrelationsData = transformCorrelationToComprehensive(correlationData, metricNames);
            end
        elseif length(firstMetric) >= 3 && ischar(firstMetric{1})
            % Single metric format [name, value, unit]
            % Transform to comprehensive format
            transformedMetricsData = transformMetricToComprehensive(metricData, metricNames);
            
            % Transform correlation data
            transformedCorrelationsData = transformCorrelationToComprehensive(correlationData, metricNames);
        else
            % Unknown format, attempt to preserve data
            transformedMetricsData = metricData;
            transformedCorrelationsData = correlationData;
        end
    else
        % Empty or invalid data, return empty results
        warning('Invalid or empty metric data format.');
    end
end

function transformedMetricsData = transformMetricToComprehensive(metricData, metricNames)
    % Transform single-metric data to comprehensive format
    transformedMetricsData = {};
    
    % If metric names not provided, use a default
    if nargin < 2 || isempty(metricNames)
        metricNames = {'Performance Metric'};
    end
    
    % Ensure metricNames is a cell array
    if ~iscell(metricNames)
        metricNames = {metricNames};
    end
    
    % Process each metric
    for i = 1:length(metricNames)
        metricName = metricNames{i};
        
        % Default values
        meanValue = 'N/A';
        stdValue = 'N/A';
        unitValue = '';
        minValue = 'N/A';
        maxValue = 'N/A';
        
        % Look for statistics in the metric data
        for j = 1:length(metricData)
            if iscell(metricData{j}) && length(metricData{j}) >= 1
                statName = metricData{j}{1};
                
                % If "Mean" statistic found for this metric
                if strcmpi(statName, 'Mean') || strcmpi(statName, 'Average')
                    if length(metricData{j}) >= 2
                        meanValue = metricData{j}{2};
                    end
                    if length(metricData{j}) >= 3
                        unitValue = metricData{j}{3};
                    end
                % If "Standard Deviation" statistic found
                elseif contains(lower(statName), 'std') || contains(lower(statName), 'deviation')
                    if length(metricData{j}) >= 2
                        stdValue = metricData{j}{2};
                    end
                % If "Minimum" statistic found
                elseif contains(lower(statName), 'min')
                    if length(metricData{j}) >= 2
                        minValue = metricData{j}{2};
                    end
                % If "Maximum" statistic found
                elseif contains(lower(statName), 'max')
                    if length(metricData{j}) >= 2
                        maxValue = metricData{j}{2};
                    end
                end
            end
        end
        
        % Add the transformed metric data
        transformedMetricsData{end+1} = {metricName, meanValue, stdValue, unitValue, minValue, maxValue};
    end
end

function transformedCorrelationsData = transformCorrelationToComprehensive(correlationData, metricNames)
    % Transform single-metric correlation data to comprehensive format
    transformedCorrelationsData = {};
    
    % If metric names not provided, use a default
    if nargin < 2 || isempty(metricNames)
        metricNames = {'Performance Metric'};
    end
    
    % Ensure metricNames is a cell array
    if ~iscell(metricNames)
        metricNames = {metricNames};
    end
    
    % Process each metric
    for i = 1:length(metricNames)
        metricName = metricNames{i};
        
        % Default correlation values
        kpCorr = 'N/A';
        kiCorr = 'N/A';
        kdCorr = 'N/A';
        
        % Look for correlations in the correlation data
        for j = 1:length(correlationData)
            if iscell(correlationData{j}) && length(correlationData{j}) >= 2
                paramName = correlationData{j}{1};
                corrValue = correlationData{j}{2};
                
                % Match parameter name to controller parameter
                if contains(lower(paramName), 'kp') || contains(lower(paramName), 'p gain')
                    kpCorr = corrValue;
                elseif contains(lower(paramName), 'ki') || contains(lower(paramName), 'i gain')
                    kiCorr = corrValue;
                elseif contains(lower(paramName), 'kd') || contains(lower(paramName), 'd gain')
                    kdCorr = corrValue;
                end
            end
        end
        
        % Add the transformed correlation data
        transformedCorrelationsData{end+1} = {metricName, kpCorr, kiCorr, kdCorr};
    end
end