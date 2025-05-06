% Script to perform batch analysis on control system parameters
% Sweeping K's numerator and denominator parameters from 0.1 to 100

% Define the plant transfer function G(s)
num_G = [-0.0004];
den_G = [1 0.1028 0.747 0.1962];
G = tf(num_G, den_G);

% Define the controller transfer function K(s)
num_K = [32 32];
den_K = [3.2 1];
K = tf(num_K, den_K);

% Define the parameter range
param_min = 0.1;
param_max = 100;
param_step = 0.5;

% Create output directory if it doesn't exist
output_dir = 'D:\Studienarbeit Regelungstechnik\BatchAnalysis';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Enable all analysis types
analysis_options = struct('stability', true, ...
                         'nyquist', true, ...
                         'bode', true, ...
                         'keyParams', true, ...
                         'margins', true, ...
                         'jump', true);

% Analyze K numerator first coefficient (index 1)
param_info = struct('type', 'K', ...
                   'coeffType', 'num', ...
                   'index', 1, ...
                   'min', param_min, ...
                   'max', param_max, ...
                   'step', param_step);

save_path = fullfile(output_dir, 'K_num_1.mat');
fprintf('Running batch analysis for K numerator coefficient 1...\n');
batchResults = runBatchAnalysisFromCommandLine(G, K, param_info, analysis_options, save_path);
fprintf('Analysis complete. Results saved to: %s\n', save_path);

% Analyze K numerator second coefficient (index 2)
param_info.index = 2;
save_path = fullfile(output_dir, 'K_num_2.mat');
fprintf('Running batch analysis for K numerator coefficient 2...\n');
batchResults = runBatchAnalysisFromCommandLine(G, K, param_info, analysis_options, save_path);
fprintf('Analysis complete. Results saved to: %s\n', save_path);

% Analyze K denominator first coefficient (index 1)
param_info.coeffType = 'den';
param_info.index = 1;
save_path = fullfile(output_dir, 'K_den_1.mat');
fprintf('Running batch analysis for K denominator coefficient 1...\n');
batchResults = runBatchAnalysisFromCommandLine(G, K, param_info, analysis_options, save_path);
fprintf('Analysis complete. Results saved to: %s\n', save_path);

% Analyze K denominator second coefficient (index 2)
param_info.index = 2;
save_path = fullfile(output_dir, 'K_den_2.mat');
fprintf('Running batch analysis for K denominator coefficient 2...\n');
batchResults = runBatchAnalysisFromCommandLine(G, K, param_info, analysis_options, save_path);
fprintf('Analysis complete. Results saved to: %s\n', save_path);

fprintf('All analyses complete!\n');