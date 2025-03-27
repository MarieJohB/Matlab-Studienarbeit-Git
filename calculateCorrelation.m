% Helper function to calculate correlation coefficient
function corr = calculateCorrelation(param, metric)
    % Remove any NaN values
    validIndices = ~isnan(param) & ~isnan(metric);
    if sum(validIndices) < 3
        corr = NaN;
        return;
    end
    
    % Calculate correlation
    corrMat = corrcoef(param(validIndices), metric(validIndices));
    if size(corrMat, 1) > 1
        corr = corrMat(1, 2);
    else
        corr = NaN;
    end
end