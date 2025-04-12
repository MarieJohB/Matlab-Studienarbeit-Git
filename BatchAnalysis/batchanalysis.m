function batchanalysis(G, K, paramType, paramCoeffType, paramIndex, paramMin, paramMax, paramStep, analysisTypes, savePath)
% BATCHANALYSIS - Command line interface to the batch analysis functionality
%
% Usage:
%   batchanalysis(G, K, paramType, paramCoeffType, paramIndex, paramMin, paramMax, paramStep, analysisTypes, savePath)
%
% Inputs:
%   G - Plant transfer function
%   K - Controller transfer function
%   paramType - 'G' or 'K' (which transfer function to modify)
%   paramCoeffType - 'num' or 'den' (which coefficient type to modify)
%   paramIndex - Index of coefficient to modify
%   paramMin - Minimum parameter value
%   paramMax - Maximum parameter value
%   paramStep - Parameter step size
%   analysisTypes - String or cell array of strings with analysis types:
%                  'stability', 'nyquist', 'bode', 'keyParams', 'jump', 'margins', 'all'
%   savePath - Path to save results
%
% Examples:
%   G = tf([1], [1 2 1]);
%   K = tf([5 1], [1 0]);
%   
%   % Sweep the proportional gain in K from 1 to 10 with step 0.5
%   batchanalysis(G, K, 'K', 'num', 1, 1, 10, 0.5, 'all', 'batch_results.mat')
%
%   % Sweep the time constant in G and analyze stability and Nyquist only
%   batchanalysis(G, K, 'G', 'den', 2, 0.1, 5, 0.1, {'stability', 'nyquist'}, 'tc_analysis.mat')

% Create parameter info structure
paramInfo = struct('type', paramType, ...
                   'coeffType', paramCoeffType, ...
                   'index', paramIndex, ...
                   'min', paramMin, ...
                   'max', paramMax, ...
                   'step', paramStep);

% Create analysis options structure
if ischar(analysisTypes)
    if strcmpi(analysisTypes, 'all')
        % All options
        analysisOptions = struct('stability', true, ...
                                'nyquist', true, ...
                                'bode', true, ...
                                'keyParams', true, ...
                                'margins', true, ...
                                'jump', true);
    else
        % Single option
        analysisOptions = struct('stability', false, ...
                                'nyquist', false, ...
                                'bode', false, ...
                                'keyParams', false, ...
                                'margins', false, ...
                                'jump', false);
        analysisOptions.(lower(analysisTypes)) = true;
    end
elseif iscell(analysisTypes)
    % Multiple options
    analysisOptions = struct('stability', false, ...
                           'nyquist', false, ...
                           'bode', false, ...
                           'keyParams', false, ...
                           'margins', false, ...
                           'jump', false);
    
    for i = 1:length(analysisTypes)
        if strcmpi(analysisTypes{i}, 'all')
            analysisOptions.stability = true;
            analysisOptions.nyquist = true;
            analysisOptions.bode = true;
            analysisOptions.keyParams = true;
            analysisOptions.margins = true;
            analysisOptions.jump = true;
            break;
        else
            analysisOptions.(lower(analysisTypes{i})) = true;
        end
    end
else
    error('analysisTypes must be a string or cell array of strings');
end

% Ensure bode is enabled if margins are requested
if analysisOptions.margins && ~analysisOptions.bode
    warning('Enabling Bode analysis because margins calculation requires it');
    analysisOptions.bode = true;
end

% Run batch analysis (no app needed for command line)
batchResults = runBatchAnalysis([], G, K, paramInfo, analysisOptions, savePath);

% Ask if user wants to visualize results
%resp = input('Do you want to visualize the results? (y/n): ', 's');
%if strcmpi(resp, 'y')
%    batchVisualization(batchResults, savePath);
%end
end