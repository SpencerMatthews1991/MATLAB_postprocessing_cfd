%% AoA Sweep Analysis for CFD Data
% This file will automatically access all angle of attack (AoA) directories
% from the file of your choosing. You will need to edit the file header and
% AoA values in "Define directories and AoA values" section. Properties are
% aribrary for Re=1000, but may be modified to your liking. The loop will
% process CFD data and calculate various performance metrics. These will be
% printed for every AoA and graphs will be generated. A comparison of raw
% CFD data versus the calculated values will be conducted. This information 
% will then be automatically saved to the file of your choosing. Besure to 
% name and date the file(s)! Cheers!!
%
% Version 2: This version imports and compares data from the Ilio (2018) 
% LBM paper. This version also generates the signal trace for the lift
% coefficicent, if you have this data. You may need to modify certain
% aspects of the post processing script to manage your data and maintain
% quality. Instructions are placed a sensitive locations to help you modify
% the inputs.
% 
% Version 3: Various quality of life updates.
% 
% Version 4: Improved metric tracking with "target" values
% 
% Version 5: More bug fixes and adjustments to increase the accuracy of
% realistic results and enhance visualization.

clc;
clear;
close all;

%% Properties of Air (standard) and CFD inputs
u = 1;          %[m/s] free stream velocity
lc = 1;         %[m] Chord length
s = 1;          %[m] span length
rho = 1.225;    %[kg/m³] density of air
nu = 0.001;     %Arbitrary Nu for analysis purposes
A = s*lc;       %[m^2] wing area
Re = (u*lc)/nu; %Reynolds number, nondimensional velocity metric

%% Define directories and AoA values
% Option 1: Specify directories manually
baseDir = './NonPorous_NACA0012';  % Update this path
aoaDirs = {'AoA0','AoA1', 'AoA2','AoA3', 'AoA4','AoA5', 'AoA6','AoA7',...
    'AoA8','AoA9', 'AoA10','AoA11','AoA12','AoA13','AoA14','AoA15','AoA16',...
    'AoA17','AoA18','AoA19','AoA20','AoA21','AoA22','AoA23','AoA24','AoA25',...
    'AoA26','AoA27','AoA28'};
aoaValues = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,...
    24,25,26,27,28];  % Corresponding AoA values in degrees

% Option 2: Automatically find directories (if they follow a pattern)
% aoaDirs = dir(fullfile(baseDir, 'AoA_*'));
% aoaDirs = {aoaDirs([aoaDirs.isdir]).name};

%% Import Ilio (2018) paper data
% Replace with alternative paper data, if you are using another paper
Il_cd = readmatrix('./ilio_2018_paper_data/Ilio_CD_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Il_cl = readmatrix('./ilio_2018_paper_data/Ilio_CL_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Il_clcd = readmatrix('./ilio_2018_paper_data/Ilio_CLCD_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Ku_cd = readmatrix('./ilio_2018_paper_data/Kurtulus_CD_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Ku_cl = readmatrix('./ilio_2018_paper_data/Kurtulus_CL_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Ku_clcd = readmatrix('./ilio_2018_paper_data/Kurtulus_CLCD_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Kh_cl = readmatrix('./ilio_2018_paper_data/Khalid_CL_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);
Li_cl = readmatrix('./ilio_2018_paper_data/Liu_CL_data.csv', 'FileType', 'text', 'NumHeaderLines', 1);

%% Import signal trace data
% Coefficient of lift signal trace data for individual runs. Replace as
% needed with your own data.
aoa4 = readmatrix('./NonPorous_NACA0012/AoA4/lift_coefficient-rfile.out', 'FileType', 'text', 'NumHeaderLines', 3);
aoa10 = readmatrix('./NonPorous_NACA0012/AoA10/lift_coefficient-rfile.out', 'FileType', 'text', 'NumHeaderLines', 3);
aoa16 = readmatrix('./NonPorous_NACA0012/AoA16/lift_coefficient-rfile.out', 'FileType', 'text', 'NumHeaderLines', 3);
aoa28 = readmatrix('./NonPorous_NACA0012/AoA28/lift_coefficient-rfile.out', 'FileType', 'text', 'NumHeaderLines', 3);

%% Plot replication from paper
% Replicates the critical plots from the Ilio (2018) paper to be used for
% later comparison with collected CFD and LBM data
figure('Position', [200 200 1200 900])

% Cl vs AoA
subplot(2,2,1)
plot(Il_cl(:,1), Il_cl(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Kh_cl(:,1),Kh_cl(:,2) , 'm-*', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Li_cl(:,1),Li_cl(:,2) , 'b-^', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Ku_cl(:,1),Ku_cl(:,2) , 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient vs AoA, Ilio(2018)')
legend('Ilio', 'Khalid', 'Liu', 'Kurtulus', 'Location', 'best')
grid on

% Cd vs AoA
subplot(2,2,3)
plot(Il_cd(:,1), Il_cd(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Ku_cd(:,1),Ku_cd(:,2) , 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Drag Coefficient (C_d)')
title('Drag Coefficient vs AoA, Ilio(2018)')
legend('Ilio', 'Kurtulus', 'Location', 'best')
grid on

% Cl/Cd vs AoA
subplot(2,2,4)
plot(Il_clcd(:,1), Il_clcd(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Ku_clcd(:,1),Ku_clcd(:,2) , 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Aerodynamic Efficiency (C_l/C_d)')
title('C_l/C_d vs AoA, Ilio(2018)')
legend('Ilio', 'Kurtulus', 'Location', 'best')
grid on

saveas(gcf,'Ilio2018_reconstruction.png','png') %Rename me! (if you want)

%% Initialize storage arrays
numCases = length(aoaValues);
Cl_calc = zeros(numCases, 1);
Cd_calc = zeros(numCases, 1);
ClCd_calc = zeros(numCases, 1);
Cl_solver = zeros(numCases, 1);
Cd_solver = zeros(numCases, 1);
ClCd_solver = zeros(numCases, 1);
Cf_solver = zeros(numCases, 1);

%% Loop through each AoA case to generate time averaged CFD data
numValues = 10000; % Number of values to time-average

for i = 1:numCases
    currentDir = fullfile(baseDir, aoaDirs{i});
    fprintf('\nProcessing AoA = %.1f degrees (%s)\n', aoaValues(i), aoaDirs{i});
    
    try
        % Load data files
        dragf = readmatrix(fullfile(currentDir, 'drag-rfile.out'), ...
            'FileType', 'text', 'NumHeaderLines', 3);
        liftf = readmatrix(fullfile(currentDir, 'lift-rfile.out'), ...
            'FileType', 'text', 'NumHeaderLines', 3);
        skinf = readmatrix(fullfile(currentDir, 'skin_friction-rfile.out'), ...
            'FileType', 'text', 'NumHeaderLines', 3);
        clcdf = readmatrix(fullfile(currentDir, 'cl-cd-rfile.out'), ...
            'FileType', 'text', 'NumHeaderLines', 3);
        clf = readmatrix(fullfile(currentDir, 'lift_coefficient-rfile.out'), ...
            'FileType', 'text', 'NumHeaderLines', 3);
        cdf = readmatrix(fullfile(currentDir, 'drag_coefficient-rfile.out'), ...
            'FileType', 'text', 'NumHeaderLines', 3);
        
        % Extract data
        dragData = dragf(:, 1:2);
        liftData = liftf(:, 1:2);
        skinData = skinf(:, 1:2);
        clcdData = clcdf(:, 1:2);
        clData = clf(:, 1:2);
        cdData = cdf(:, 1:2);
        
        % Time averaging (last numValues)
        if size(dragData, 1) > numValues
            dragData = dragData(end-numValues+1:end, :);
        end
        if size(liftData, 1) > numValues
            liftData = liftData(end-numValues+1:end, :);
        end
        if size(skinData, 1) > numValues
            skinData = skinData(end-numValues+1:end, :);
        end
        if size(clcdData, 1) > numValues
            clcdData = clcdData(end-numValues+1:end, :);
        end
        if size(clData, 1) > numValues
            clData = clData(end-numValues+1:end, :);
        end
        if size(cdData, 1) > numValues
            cdData = cdData(end-numValues+1:end, :);
        end
        
        % Calculate time averages
        time_avgDragF = mean(dragData(:, 2));
        time_avgLiftF = mean(liftData(:, 2));
        time_avgSkinF = mean(skinData(:, 2));
        time_avgClCdF = mean(clcdData(:, 2));
        time_avgClF = mean(clData(:, 2));
        time_avgCdF = mean(cdData(:, 2));
        
        % Calculate coefficients
        Cd_calc(i) = (2 * time_avgDragF) / (rho * u^2 * A);
        Cl_calc(i) = (2 * time_avgLiftF) / (rho * u^2 * A);
        %ClCd_calc(i) = Cl_calc(i) / Cd_calc(i);
        Cl_solver(i) = time_avgClF;
        Cd_solver(i) = time_avgCdF;
        ClCd_solver(i) = time_avgClCdF;
        Cf_solver(i) = time_avgSkinF;
        ClCd_calc(i) = Cl_solver(i) / Cd_solver(i);
        
        % Display results
        fprintf('  Cl (calc): %.4f, Cl (CFD): %.4f\n', Cl_calc(i), Cl_solver(i));
        fprintf('  Cd (calc): %.4f, Cd (CFD): %.4f\n', Cd_calc(i), Cd_solver(i));
        fprintf('  Cl/Cd (calc): %.4f, Cl/Cd (CFD): %.4f\n', ClCd_calc(i), ClCd_solver(i));
        
    catch ME
        warning('Error processing AoA = %.1f: %s', aoaValues(i), ME.message);
        % Fill with NaN if processing fails
        Cl_calc(i) = NaN;
        Cd_calc(i) = NaN;
        ClCd_calc(i) = NaN;
        Cl_solver(i) = NaN;
        Cd_solver(i) = NaN;
        ClCd_solver(i) = NaN;
    end
end

%% Generate AoA sweep plots
figure('Position', [200 200 1200 900])

% Cl vs AoA
subplot(2,2,1)
plot(Ku_cl(:,1), Ku_cl(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Il_cl(:,1), Il_cl(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, Cl_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, Cl_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)')
ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient vs AoA')
legend('Target: Kurtulus', 'Target: Ilio', 'Calculated', 'CFD Solver', 'Location', 'best')
grid on
hold off

% Cd vs AoA
subplot(2,2,2)
plot(Ku_cd(:,1), Ku_cd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Il_cd(:,1), Il_cd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, Cd_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, Cd_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)')
ylabel('Drag Coefficient (C_d)')
title('Drag Coefficient vs AoA')
legend('Target: Kurtulus', 'Target: Ilio', 'Calculated', 'CFD Solver', 'Location', 'best')
grid on

% Cl/Cd vs AoA
subplot(2,2,3)
plot(Ku_clcd(:,1), Ku_clcd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Il_clcd(:,1), Il_clcd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, ClCd_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, ClCd_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)')
ylabel('Aerodynamic Efficiency (C_l/C_d)')
title('C_l/C_d vs AoA')
legend('Target: Kurtulus', 'Target: Ilio', 'Calculated', 'CFD Solver', 'Location', 'best')
grid on

% Cl vs Cd 
subplot(2,2,4)
plot(Cd_calc, Cl_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Cd_solver, Cl_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Drag Coefficient (C_d)')
ylabel('Lift Coefficient (C_l)')
title('C_l vs C_d')
legend('Calculated', 'CFD Solver', 'Location', 'best')
grid on

saveas(gcf,'20260214_SM_dataprocess_out.png','png') %Rename me!

% Find optimal AoA
[maxClCd, maxIdx] = max(ClCd_solver);
fprintf('\nOptimal AoA: %.1f degrees with Cl/Cd = %.4f\n', aoaValues(maxIdx), maxClCd);

%% Save results
resultsTable = table(aoaValues', Cl_calc, Cd_calc, ClCd_calc, ...
    Cl_solver, Cd_solver, ClCd_solver, ...
    'VariableNames', {'AoA_deg', 'Cl_calc', 'Cd_calc', 'ClCd_calc', ...
    'Cl_solver', 'Cd_solver', 'ClCd_solver'});
writetable(resultsTable, '20260214_AoA_sweep_results.csv'); %Rename me!
fprintf('\nResults saved to AoA_sweep_results.csv\n');

%% Plot Comparison of CFD data with previous research data
figure('Position', [200 200 1200 900])

% Cl vs AoA
subplot(2,2,1)
plot(Il_cl(:,1), Il_cl(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Kh_cl(:,1),Kh_cl(:,2) , 'm-*', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Li_cl(:,1),Li_cl(:,2) , 'b-^', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Ku_cl(:,1),Ku_cl(:,2) , 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, Cl_solver, 'c-v', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient (C_l) vs AoA')
legend('Ilio', 'Khalid', 'Liu', 'Kurtulus', 'Matthews', 'Location', 'best')
grid on

% CL solver signal trace
subplot(2,2,2)
plot(aoa4(:,3), aoa4(:,2), 'b-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoa10(:,3), aoa10(:,2), 'm-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoa16(:,3), aoa16(:,2), 'r-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoa28(:,3), aoa28(:,2), 'g-', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Time Step Iteration')
ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient (C_l) CFD Signal Trace')
legend('AoA = 4', 'AoA = 10', 'AoA = 16', 'AoA = 28', 'Location', 'best')
ylim ([0 2.5])
grid on

% Cd vs AoA
subplot(2,2,3)
plot(Il_cd(:,1), Il_cd(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Ku_cd(:,1),Ku_cd(:,2) , 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, Cd_solver, 'c-v', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Drag Coefficient (C_d)')
title('Drag Coefficient (C_d) vs AoA')
legend('Ilio', 'Kurtulus', 'Matthews', 'Location', 'best')
grid on

% % Cl/Cd vs AoA
subplot(2,2,4)
plot(Il_clcd(:,1), Il_clcd(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(Ku_clcd(:,1),Ku_clcd(:,2) , 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
hold on
plot(aoaValues, ClCd_calc, 'c-v', 'LineWidth', 2, 'MarkerSize', 8)
hold on
xlabel('Angle of Attack (degrees)')
ylabel('Aerodynamic Efficiency (C_l/C_d)')
title('Aerodynamic Efficiency (C_l/C_d) vs AoA')
legend('Ilio', 'Kurtulus', 'Matthews', 'Location', 'best')
grid on

saveas(gcf,'20260214_SM_dataprocess_comparison.png','png') %Rename me!

%% Calculate overall difference between Matthews and Ilio data
% Interpolate Ilio data to match AoA values
Il_cl_interp = interp1(Il_cl(:,1), Il_cl(:,2), aoaValues, 'linear', 'extrap');
Il_cd_interp = interp1(Il_cd(:,1), Il_cd(:,2), aoaValues, 'linear', 'extrap');
Il_clcd_interp = interp1(Il_clcd(:,1), Il_clcd(:,2), aoaValues, 'linear', 'extrap');

% Calculate point-wise percent differences
Cl_percentDiff = ((Cl_solver' - Il_cl_interp) ./ Il_cl_interp) * 100;
Cl_percentDiff(1) = NaN;    % Omit first value close to zero
Cd_percentDiff = ((Cd_solver' - Il_cd_interp) ./ Il_cd_interp) * 100;
ClCd_percentDiff = ((ClCd_calc' - Il_clcd_interp) ./ Il_clcd_interp) * 100;
ClCd_percentDiff(1) = NaN;  % Omit first value divided by NaN

% Calculate mean absolute percent difference (overall metric)
Cl_MAPD = mean(abs(Cl_percentDiff), 'omitnan');
Cd_MAPD = mean(abs(Cd_percentDiff), 'omitnan');
ClCd_MAPD = mean(abs(ClCd_percentDiff), 'omitnan');

% Calculate RMS percent difference (alternative metric)
Cl_RMSPD = sqrt(mean(Cl_percentDiff.^2, 'omitnan'));
Cd_RMSPD = sqrt(mean(Cd_percentDiff.^2, 'omitnan'));
ClCd_RMSPD = sqrt(mean(ClCd_percentDiff.^2, 'omitnan'));

% Display results
fprintf('\n========== Overall Comparison: Matthews vs Ilio ==========\n');
fprintf('Lift Coefficient (Cl):\n');
fprintf('  Mean Absolute Percent Difference: %.2f%%\n', Cl_MAPD);
fprintf('  RMS Percent Difference: %.2f%%\n', Cl_RMSPD);
fprintf('\nDrag Coefficient (Cd):\n');
fprintf('  Mean Absolute Percent Difference: %.2f%%\n', Cd_MAPD);
fprintf('  RMS Percent Difference: %.2f%%\n', Cd_RMSPD);
fprintf('\nAerodynamic Efficiency (Cl/Cd):\n');
fprintf('  Mean Absolute Percent Difference: %.2f%%\n', ClCd_MAPD);
fprintf('  RMS Percent Difference: %.2f%%\n', ClCd_RMSPD);
fprintf('==========================================================\n');