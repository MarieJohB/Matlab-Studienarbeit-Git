function createSummaryTab(tab, batchResults)
% Create a summary tab with HTML tables for better visual display and
% integrated rendered transfer functions - without panels

% Extract parameter info
paramInfo = batchResults.paramInfo;
paramValues = batchResults.paramValues;

% Get transfer function data for rendering
G = batchResults.baseG;
K = batchResults.baseK;
[G_num, G_den] = tfdata(G, 'v');
[K_num, K_den] = tfdata(K, 'v');

% Create transfer function HTML
G_html = polyToHTMLString(G_num);
G_den_html = polyToHTMLString(G_den);
K_html = polyToHTMLString(K_num);
K_den_html = polyToHTMLString(K_den);

% Create batch analysis info as HTML table with integrated transfer functions
htmlContent = ['<html><head><style>', ...
    'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; }', ...
    'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 14px; }', ...
    'td { border: 1px solid #DDDDDD; padding: 8px; font-size: 13px; }', ...
    'tr:nth-child(even) { background-color: #F2F2F2; }', ...
    'tr:hover { background-color: #E8F4F8; }', ...
    '.header { background-color: #5B9BD5; color: white; font-weight: bold; }', ...
    '.fraction { display: inline-block; vertical-align: middle; text-align: center; }', ...
    '.fraction .num { border-bottom: 1.5px solid black; padding: 0px 3px; }', ...
    '.fraction .den { padding: 0px 3px; }', ...
    '</style></head><body>', ...
    '<table>', ...
    '<tr class="header"><th colspan="2">Batch Analysis Information</th></tr>'];

% Parameter info
paramStr = sprintf('%s: %s [%d]', paramInfo.type, paramInfo.coeffType, paramInfo.index);
minVal = paramInfo.min;
maxVal = paramInfo.max;
stepSize = paramInfo.step;
numValues = length(batchResults.paramValues);

% Add parameter info rows to the table
htmlContent = [htmlContent, ...
    '<tr><td>Parameter Swept</td><td>' paramStr '</td></tr>', ...
    '<tr><td>Range</td><td>' sprintf('%.4f to %.4f with step size %.4f (%d values)', ...
        minVal, maxVal, stepSize, numValues) '</td></tr>'];

% Add G(s) with rendered equation
htmlContent = [htmlContent, ...
    '<tr><td>Base Plant</td><td>G(s) = <div class="fraction"><div class="num">' G_html '</div><div class="den">' G_den_html '</div></div></td></tr>'];

% Add K(s) with rendered equation
htmlContent = [htmlContent, ...
    '<tr><td>Base Controller</td><td>K(s) = <div class="fraction"><div class="num">' K_html '</div><div class="den">' K_den_html '</div></div></td></tr>', ...
    '</table></body></html>'];

% Create HTML component for batch info - directly on tab
batchInfoHtml = uihtml(tab, 'HTMLSource', htmlContent, 'Position', [50 735 1800 220]);

% Add title label for batch info
uilabel(tab, 'Position', [50 970 400 20], 'Text', 'Batch Analysis Information', ...
    'FontWeight', 'bold', 'FontSize', 16);

% Create key findings HTML table if stability/margins info is available
if isfield(batchResults, 'stability') || isfield(batchResults, 'margins')
    findingsHtml = createKeyFindingsHTML(batchResults);
    uihtml(tab, 'HTMLSource', findingsHtml, 'Position', [50 450 1800 290]);
end

% Plot parameter vs stability if available (larger size for fullscreen)
if isfield(batchResults, 'stability')
    stabilityAxes = uiaxes(tab, 'Position', [50 50 1800 380]);
    plot(stabilityAxes, batchResults.paramValues, double(batchResults.stability), ...
        'LineWidth', 3, 'Marker', '.', 'MarkerSize', 20, 'Color', [0.3 0.6 0.9]);
    title(stabilityAxes, 'Parameter vs. Stability', 'FontSize', 16, 'FontWeight', 'bold');
    xlabel(stabilityAxes, 'Parameter Value', 'FontSize', 14);
    ylabel(stabilityAxes, 'Stability (0=Unstable, 1=Stable)', 'FontSize', 14);
    ylim(stabilityAxes, [-0.1 1.1]);
    grid(stabilityAxes, 'on');
    
    % Add stability transitions markers if any
    transitions = findStabilityTransitions(batchResults);
    if ~isempty(transitions)
        hold(stabilityAxes, 'on');
        for t = transitions
            transValue = (paramValues(t) + paramValues(t+1))/2;
            xline(stabilityAxes, transValue, 'r--', 'LineWidth', 2);
            text(stabilityAxes, transValue, 0.5, 'Stability Transition', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontWeight', 'bold', 'Color', 'r', 'FontSize', 12);
        end
        hold(stabilityAxes, 'off');
    end
end
end

function htmlContent = createKeyFindingsHTML(batchResults)
% Create HTML content for key findings

htmlContent = ['<html><head><style>', ...
    'table { border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; }', ...
    'th { background-color: #4472C4; color: white; padding: 8px; text-align: center; font-size: 14px; }', ...
    'td { border: 1px solid #DDDDDD; padding: 8px; font-size: 13px; }', ...
    'tr:nth-child(even) { background-color: #F2F2F2; }', ...
    'tr:hover { background-color: #E8F4F8; }', ...
    '.header { background-color: #5B9BD5; color: white; font-weight: bold; }', ...
    '.good { color: green; font-weight: bold; }', ...
    '.bad { color: red; font-weight: bold; }', ...
    '.warning { color: orange; font-weight: bold; }', ...
    '</style></head><body>', ...
    '<table>', ...
    '<tr class="header"><th colspan="2">Key Findings</th></tr>'];

% Add stability findings if available
if isfield(batchResults, 'stability')
    % Calculate stability percentage
    stabilityPercent = sum(batchResults.stability) / length(batchResults.stability) * 100;
    
    % Calculate stability transitions
    transitions = findStabilityTransitions(batchResults);
    
    if isempty(transitions)
        if all(batchResults.stability)
            htmlContent = [htmlContent, ...
                '<tr><td>Stability</td><td class="good">System remains stable across the entire parameter range.</td></tr>'];
        elseif ~any(batchResults.stability)
            htmlContent = [htmlContent, ...
                '<tr><td>Stability</td><td class="bad">System is unstable across the entire parameter range.</td></tr>'];
        end
    else
        if stabilityPercent >= 75
            stabilityClass = 'good';
        elseif stabilityPercent >= 25
            stabilityClass = 'warning';
        else
            stabilityClass = 'bad';
        end
        
        htmlContent = [htmlContent, ...
            '<tr><td>Stability</td><td class="' stabilityClass '">' ...
            sprintf('%.1f%% of parameter values yield a stable system. ', stabilityPercent) ...
            sprintf('Found %d stability transition(s).', length(transitions)) ...
            '</td></tr>'];
            
        % Add first few transitions
        for i = 1:min(3, length(transitions))
            idx = transitions(i);
            transValue = (batchResults.paramValues(idx) + batchResults.paramValues(idx+1))/2;
            if batchResults.stability(idx)
                transType = 'Stable → Unstable';
            else
                transType = 'Unstable → Stable';
            end
            
            htmlContent = [htmlContent, ...
                '<tr><td>Transition ' num2str(i) '</td><td>' ...
                sprintf('At parameter = %.4f: %s', transValue, transType) ...
                '</td></tr>'];
        end
    end
end

% Add margins findings if available
if isfield(batchResults, 'margins')
    % Initialize arrays for margins
    numPoints = length(batchResults.paramValues);
    gainMargin_dB = zeros(1, numPoints);
    phaseMargin_deg = zeros(1, numPoints);
    
    % Extract data for each parameter value
    for i = 1:numPoints
        if ~isempty(batchResults.margins{i})
            margins = batchResults.margins{i};
            
            % Convert gain margin to dB
            if ~isnan(margins.gainMargin)
                gainMargin_dB(i) = 20*log10(margins.gainMargin);
            else
                gainMargin_dB(i) = NaN;
            end
            
            phaseMargin_deg(i) = margins.phaseMargin;
        else
            gainMargin_dB(i) = NaN;
            phaseMargin_deg(i) = NaN;
        end
    end
    
    % Count values meeting criteria
    minGainMargin_dB = 6.0;  % dB
    minPhaseMargin_deg = 30.0;  % degrees
    
    goodGmCount = sum(gainMargin_dB >= minGainMargin_dB & ~isnan(gainMargin_dB));
    goodPmCount = sum(phaseMargin_deg >= minPhaseMargin_deg & ~isnan(phaseMargin_deg));
    bothGoodCount = sum((gainMargin_dB >= minGainMargin_dB & phaseMargin_deg >= minPhaseMargin_deg) & ...
        ~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));
    validCount = sum(~isnan(gainMargin_dB) & ~isnan(phaseMargin_deg));
    
    if validCount > 0
        bothGoodPercent = bothGoodCount / validCount * 100;
        
        if bothGoodPercent >= 75
            marginClass = 'good';
        elseif bothGoodPercent >= 25
            marginClass = 'warning';
        else
            marginClass = 'bad';
        end
        
        htmlContent = [htmlContent, ...
            '<tr><td>Stability Margins</td><td class="' marginClass '">' ...
            sprintf('%.1f%% of parameter values have adequate gain and phase margins.', bothGoodPercent) ...
            '</td></tr>'];
            
        % Find best parameter for margins
        if bothGoodCount > 0
            % Find parameter value with best compromise (highest combined normalized margins)
            normGm = gainMargin_dB / minGainMargin_dB;
            normPm = phaseMargin_deg / minPhaseMargin_deg;
            combinedNorm = normGm + normPm;
            [~, bestIdx] = max(combinedNorm);
            
            htmlContent = [htmlContent, ...
                '<tr><td>Best Parameter</td><td>' ...
                sprintf('Parameter = %.4f provides GM = %.2f dB, PM = %.2f°', ...
                batchResults.paramValues(bestIdx), gainMargin_dB(bestIdx), phaseMargin_deg(bestIdx)) ...
                '</td></tr>'];
        end
    end
end

htmlContent = [htmlContent, '</table></body></html>'];
end

function polyStr = polyToHTMLString(coeff)
% Convert a polynomial coefficient vector to an HTML string
% with proper formatting for display

% Remove any leading zeros to get the correct degree
idx = find(coeff ~= 0, 1, 'first');
if isempty(idx)
    polyStr = '0';
    return;
end
coeff = coeff(idx:end);

deg = length(coeff) - 1;  % Determine the degree of the polynomial
terms = {};               % Cell array to store individual terms

% Iterate over all coefficients
for i = 1:length(coeff)
    coef = coeff(i);
    exp = deg - (i - 1);  % Determine the exponent for the current term
    
    % Skip zero coefficients as they don't affect the expression
    if coef == 0
        continue;
    end
    
    % Handle signs: add " + " if not the first term;
    % For negative coefficients, add " - " and use the absolute value.
    if coef > 0 && ~isempty(terms)
        term = ' + ';
    elseif coef < 0
        term = ' - ';
        coef = abs(coef);  % Use absolute value for display
    else
        term = '';
    end
    
    % Show the coefficient unless it's 1 and not the constant term (exp==0)
    if coef ~= 1 || exp == 0
        term = strcat(term, num2str(coef));
    end
    
    % Add "s" and superscript for exponents greater than 0.
    if exp > 1
        term = strcat(term, 's<sup>', num2str(exp), '</sup>');
    elseif exp == 1
        term = strcat(term, 's');
    end
    
    % Add the formatted term to the cell array.
    terms{end + 1} = term;
end

% Join all terms into a single string
polyStr = strjoin(terms, '');

% If all coefficients are zero, return "0"
if isempty(polyStr)
    polyStr = '0';
end
end