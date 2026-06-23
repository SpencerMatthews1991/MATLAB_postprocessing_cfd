clc;
clear;
close all;

%% Configuration
numValues = 10000; % Trailing rows to average (uses all rows if fewer exist)

%% Locate files relative to this script's folder
scriptDir = fileparts(mfilename('.'));
cdPath = fullfile(scriptDir, 'cd_file.csv');
clPath = fullfile(scriptDir, 'cl_file.csv');
cyPath = fullfile(scriptDir, 'cy_file.csv');

% Verify the files actually exist before trying to read them
if ~isfile(cdPath)
    error('Cannot find cd_file.csv. Expected at: %s', cdPath);
end
if ~isfile(clPath)
    error('Cannot find cl_file.csv. Expected at: %s', clPath);
end
if ~isfile(clPath)
    error('Cannot find cy_file.csv. Expected at: %s', cyPath);
end

%% Load data files
cd_data = readmatrix(cdPath, 'FileType', 'text', 'NumHeaderLines', 1);
cl_data = readmatrix(clPath, 'FileType', 'text', 'NumHeaderLines', 1);
cy_data = readmatrix(cyPath, 'FileType', 'text', 'NumHeaderLines', 1);

% col 1 = Iteration, col 2 = Cd or Cl value

%% Trim to last numValues rows (if data is longer than numValues)
if size(cd_data, 1) > numValues
    cd_data = cd_data(end-numValues+1:end, :);
end
if size(cl_data, 1) > numValues
    cl_data = cl_data(end-numValues+1:end, :);
end
if size(cy_data, 1) > numValues
    cy_data = cy_data(end-numValues+1:end, :);
end

%% Compute time averages
Cd_avg = mean(cd_data(:, 2));
Cl_avg = mean(cl_data(:, 2));
Cy_avg = mean(cy_data(:, 2));

%% Display results
fprintf('--- Time-Averaged Results ---\n');
fprintf('Rows used for Cd average : %d\n', size(cd_data, 1));
fprintf('Rows used for Cl average : %d\n', size(cl_data, 1));
fprintf('Rows used for Cy average : %d\n', size(cy_data, 1));
fprintf('Time-averaged Cd         : %.6f\n', Cd_avg);
fprintf('Time-averaged Cl         : %.6f\n', Cl_avg);
fprintf('Time-averaged Cy         : %.6f\n', Cy_avg);