%% AoA Sweep Analysis for CFD Data
% This file will automatically access all angle of attack (AoA) directories
% from the file of your choosing. You will need to edit the file header and
% AoA values in "Define directories and AoA values" section. Properties are
% arbitrary for Re=1000, but may be modified to your liking. The loop will
% process CFD data and calculate various performance metrics. These will be
% printed for every AoA and graphs will be generated. A comparison of raw
% CFD data versus the calculated values will be conducted. This information
% will then be automatically saved to the file of your choosing. Be sure to
% name and date the file(s)! Cheers!!
%
% Version 2: Imports/compares Ilio (2018) LBM data; generates Cl signal trace.
% Version 3: Various quality of life updates.
% Version 4: Improved metric tracking with "target" values.
% Version 5: Bug fixes and visualization accuracy improvements.
% Version 6: Reynolds-stress and wake analysis for the full AoA sweep.
%
% Version 7 (robustness overhaul):
%   - Files are now located by CONTENT, not exact name. Each signal is matched
%     against a normalized key (lowercase, separators and the 'rfile' token
%     stripped), so 'drag-rfile.out', 'drag_rfile.out', 'Drag.out' all resolve
%     to the same signal, while 'drag' stays distinct from 'drag_coefficient'.
%   - The value column is found by READING THE FLUENT HEADER NAMES and auto-
%     detecting the number of header lines, instead of assuming columns 1:2.
%     Works for 2-column [iter, value] and 3-column [iter, value, flow-time].
%   - Each signal loads independently: one missing/misnamed file no longer
%     NaNs out the entire AoA case.
%   - Fixed bug: wake total pressure previously read from the static file.
%   See local functions at the BOTTOM of this file for the implementation.

clc;
clear;
close all;

%% Properties of Air (Modified) and CFD inputs
u = 1;          %[m/s] free stream velocity
lc = 1;         %[m] Chord length
s = 1;          %[m] span length
rho = 1;        %[kg/m³] density of air
nu = 0.001;     %Arbitrary Nu for analysis purposes
A = s*lc;       %[m^2] wing area
Re = (u*lc)/nu; %Reynolds number, nondimensional velocity metric

%% Define directories and AoA values
% Option 1: Specify directories manually
baseDir = './NonPorous_NACA0012';  % Update this path
aoaDirs = {'AoA0','AoA1', 'AoA2','AoA3', 'AoA4','AoA5', 'AoA6','AoA7',...
    'AoA8','AoA9', 'AoA10','AoA11','AoA12','AoA13','AoA14','AoA15','AoA16',...
    'AoA17','AoA18','AoA19','AoA20'};
% ,'AoA21','AoA22','AoA23','AoA23',...
%     'AoA25','AoA26','AoA27','AoA28'};
aoaValues = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    % ...
    % 21,22,23,24,25,26,27,28];

% Option 2: Automatically find directories (if they follow a pattern)
% aoaDirs = dir(fullfile(baseDir, 'AoA*'));
% aoaDirs = {aoaDirs([aoaDirs.isdir]).name};

%% ----------------------------------------------------------------------- %%
%% SIGNAL MANIFEST  (the only thing you edit to add/remove/rename signals)  %%
%% Each row:  {storage field name , {acceptable name keys, priority order}} %%
%% Name keys are matched after normalization, so you do NOT need to match   %%
%% punctuation, case, or the 'rfile' suffix -- 'uu_wake_points-rfile.out'   %%
%% and 'UUWakePoints.out' both satisfy the key 'uuwakepoints'.              %%
%% ----------------------------------------------------------------------- %%
signalManifest = {
    'drag',    {'drag','dragforce'}
    'lift',    {'lift','liftforce'}
    'skin',    {'skinfriction','cf','skinfrictioncoefficient'}
    'clcd',    {'clcd','clcdratio','liftdragratio'}
    'cl',      {'liftcoefficient','cl'}
    'cd',      {'dragcoefficient','cd'}
    'uu',      {'uu','uureynoldsstress'}
    'vv',      {'vv','vvreynoldsstress'}
    'uv',      {'uv','uvreynoldsstress'}
    'uuwake',  {'uuwakepoints','uuwake'}
    'vvwake',  {'vvwakepoints','vvwake'}
    'uvwake',  {'uvwakepoints','uvwake'}
    'uwake',   {'velocityxwakepoints','velocityxwake','uwakepoints'}
    'vwake',   {'velocityywakepoints','velocityywake','vwakepoints'}
    'vmag',    {'velocitymagwakepoints','velocitymagnitudewakepoints','vmagwakepoints'}
    'static',  {'staticpressurewakepoints','staticpressurewake'}
    'total',   {'totalpressurewakepoints','totalpressurewake'}
};

%% Import Ilio (2018) paper data
% Replace with alternative paper data, if you are using another paper
Il_cd   = readmatrix('./ilio_2018_paper_data/Ilio_CD_data.csv',     'FileType','text','NumHeaderLines',1);
Il_cl   = readmatrix('./ilio_2018_paper_data/Ilio_CL_data.csv',     'FileType','text','NumHeaderLines',1);
Il_clcd = readmatrix('./ilio_2018_paper_data/Ilio_CLCD_data.csv',   'FileType','text','NumHeaderLines',1);
Ku_cd   = readmatrix('./ilio_2018_paper_data/Kurtulus_CD_data.csv', 'FileType','text','NumHeaderLines',1);
Ku_cl   = readmatrix('./ilio_2018_paper_data/Kurtulus_CL_data.csv', 'FileType','text','NumHeaderLines',1);
Ku_clcd = readmatrix('./ilio_2018_paper_data/Kurtulus_CLCD_data.csv','FileType','text','NumHeaderLines',1);
Kh_cl   = readmatrix('./ilio_2018_paper_data/Khalid_CL_data.csv',   'FileType','text','NumHeaderLines',1);
Li_cl   = readmatrix('./ilio_2018_paper_data/Liu_CL_data.csv',      'FileType','text','NumHeaderLines',1);

%% Import signal trace data (robust)
% Cl signal traces for individual runs. Directories and colours are listed in
% traceConfig; any directory that is missing or has no Cl file is simply skipped.
traceConfig = {'AoA4','b-'; 'AoA10','m-'; 'AoA16','r-'; 'AoA28','g-'};
traces = struct();
for t = 1:size(traceConfig,1)
    d = traceConfig{t,1};
    dpath = fullfile(baseDir, d);
    tr = [];
    if isfolder(dpath)
        fmap = buildFileMap(dpath);
        fp = findReport(fmap, {'liftcoefficient','cl'});
        if ~isempty(fp)
            tr = readReport(fp);
        end
    end
    traces.(d) = tr;   % [] if unavailable, else struct with .x and .y
end

%% Plot replication from paper
% Replicates the critical plots from the Ilio (2018) paper for later comparison
figure('Position', [200 200 1200 900])

subplot(2,2,1)
plot(Il_cl(:,1), Il_cl(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Kh_cl(:,1), Kh_cl(:,2), 'm-*', 'LineWidth', 2, 'MarkerSize', 8)
plot(Li_cl(:,1), Li_cl(:,2), 'b-^', 'LineWidth', 2, 'MarkerSize', 8)
plot(Ku_cl(:,1), Ku_cl(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient vs AoA, Ilio(2018)')
legend('Ilio','Khalid','Liu','Kurtulus','Location','best'); grid on

subplot(2,2,3)
plot(Il_cd(:,1), Il_cd(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Ku_cd(:,1), Ku_cd(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Drag Coefficient (C_d)')
title('Drag Coefficient vs AoA, Ilio(2018)')
legend('Ilio','Kurtulus','Location','best'); grid on

subplot(2,2,4)
plot(Il_clcd(:,1), Il_clcd(:,2), 'k-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Ku_clcd(:,1), Ku_clcd(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Aerodynamic Efficiency (C_l/C_d)')
title('C_l/C_d vs AoA, Ilio(2018)')
legend('Ilio','Kurtulus','Location','best'); grid on

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
uu_calc = zeros(numCases,1);
uuwake_calc = zeros(numCases,1);
vv_calc = zeros(numCases,1);
vvwake_calc = zeros(numCases,1);
uv_calc = zeros(numCases,1);
uvwake_calc = zeros(numCases,1);
StaticPressure_wake = zeros(numCases,1);
TotalPressure_wake = zeros(numCases,1);
VelocityMag_wake = zeros(numCases,1);
VelocityX_wake = zeros(numCases,1);
VelocityY_wake = zeros(numCases,1);

%% Loop through each AoA case to generate time averaged CFD data
numValues = 10000; % Number of trailing values to time-average

for i = 1:numCases
    currentDir = fullfile(baseDir, aoaDirs{i});
    fprintf('\nProcessing AoA = %.1f degrees (%s)\n', aoaValues(i), aoaDirs{i});

    % Build a normalized name -> path map of every .out file in this directory.
    if ~isfolder(currentDir)
        warning('Directory not found: %s -- filling AoA %.1f with NaN.', ...
            currentDir, aoaValues(i));
        fmap = containers.Map('KeyType','char','ValueType','char');
    else
        fmap = buildFileMap(currentDir);
    end

    % Load and time-average every signal in the manifest, independently.
    avg = struct();
    for sIdx = 1:size(signalManifest,1)
        field = signalManifest{sIdx,1};
        keys  = signalManifest{sIdx,2};
        fp = findReport(fmap, keys);
        if isempty(fp)
            warning('  [%s] no file matching {%s} -- set to NaN.', ...
                field, strjoin(keys, ', '));
            avg.(field) = NaN;
            continue;
        end
        try
            R = readReport(fp);              % name-based column detection
            avg.(field) = tailAvg(R.y, numValues);
        catch ME
            warning('  [%s] failed to read %s: %s', field, fp, ME.message);
            avg.(field) = NaN;
        end
    end

    % Calculate coefficients (identical formulas to previous versions)
    Cd_calc(i) = (2 * avg.drag) / (rho * u^2 * A);
    Cl_calc(i) = (2 * avg.lift) / (rho * u^2 * A);
    Cl_solver(i)   = avg.cl;
    Cd_solver(i)   = avg.cd;
    ClCd_solver(i) = avg.clcd;
    Cf_solver(i)   = avg.skin;
    ClCd_calc(i)   = Cl_solver(i) / Cd_solver(i);
    uu_calc(i) = avg.uu;
    vv_calc(i) = avg.vv;
    uv_calc(i) = avg.uv;
    uuwake_calc(i) = avg.uuwake;
    vvwake_calc(i) = avg.vvwake;
    uvwake_calc(i) = avg.uvwake;
    VelocityX_wake(i)      = avg.uwake;
    VelocityY_wake(i)      = avg.vwake;
    VelocityMag_wake(i)    = avg.vmag;
    StaticPressure_wake(i) = avg.static;
    TotalPressure_wake(i)  = avg.total;   % FIXED: own file, not the static file

    % Display results
    fprintf('  Cl (calc): %.4f, Cl (CFD): %.4f\n', Cl_calc(i), Cl_solver(i));
    fprintf('  Cd (calc): %.4f, Cd (CFD): %.4f\n', Cd_calc(i), Cd_solver(i));
    fprintf('  Cl/Cd (calc): %.4f, Cl/Cd (CFD): %.4f\n', ClCd_calc(i), ClCd_solver(i));
end

%% Generate AoA sweep plots
figure('Position', [200 200 1200 900])

subplot(2,2,1)
plot(Ku_cl(:,1), Ku_cl(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Il_cl(:,1), Il_cl(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, Cl_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, Cl_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient vs AoA')
legend('Target: Kurtulus','Target: Ilio','Calculated','CFD Solver','Location','best')
grid on; hold off

subplot(2,2,2)
plot(Ku_cd(:,1), Ku_cd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Il_cd(:,1), Il_cd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, Cd_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, Cd_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Drag Coefficient (C_d)')
title('Drag Coefficient vs AoA')
legend('Target: Kurtulus','Target: Ilio','Calculated','CFD Solver','Location','best')
grid on

subplot(2,2,3)
plot(Ku_clcd(:,1), Ku_clcd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Il_clcd(:,1), Il_clcd(:,2), '-x', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, ClCd_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, ClCd_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Aerodynamic Efficiency (C_l/C_d)')
title('C_l/C_d vs AoA')
legend('Target: Kurtulus','Target: Ilio','Calculated','CFD Solver','Location','best')
grid on

subplot(2,2,4)
plot(Cd_calc, Cl_calc, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Cd_solver, Cl_solver, 'b--s', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Drag Coefficient (C_d)'); ylabel('Lift Coefficient (C_l)')
title('C_l vs C_d')
legend('Calculated','CFD Solver','Location','best'); grid on

saveas(gcf,'20260624_SM_dataprocess_out.png','png') %Rename me!

% Find optimal AoA
[maxClCd, maxIdx] = max(ClCd_calc);
fprintf('\nOptimal AoA: %.1f degrees with Cl/Cd = %.4f\n', aoaValues(maxIdx), maxClCd);

%% Save results
resultsTable = table(aoaValues', Cl_calc, Cd_calc, ClCd_calc, ...
    Cl_solver, Cd_solver, ClCd_solver, ...
    'VariableNames', {'AoA_deg','Cl_calc','Cd_calc','ClCd_calc', ...
    'Cl_solver','Cd_solver','ClCd_solver'});
writetable(resultsTable, '20260624_AoA_sweep_results.csv'); %Rename me!
fprintf('\nResults saved to 20260624_AoA_sweep_results.csv\n');

%% Plot Comparison of CFD data with previous research data
figure('Position', [200 200 1200 900])

subplot(2,2,1)
plot(Il_cl(:,1), Il_cl(:,2), 'w-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Kh_cl(:,1), Kh_cl(:,2), 'm-*', 'LineWidth', 2, 'MarkerSize', 8)
plot(Li_cl(:,1), Li_cl(:,2), 'b-^', 'LineWidth', 2, 'MarkerSize', 8)
plot(Ku_cl(:,1), Ku_cl(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, Cl_solver, 'c-v', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient (C_l) vs AoA')
legend('Ilio','Khalid','Liu','Kurtulus','Matthews','Location','best'); grid on

% CL solver signal trace (robust: plots whichever traces were found)
subplot(2,2,2)
hold on
legEntries = {};
for t = 1:size(traceConfig,1)
    d  = traceConfig{t,1};
    tr = traces.(d);
    if ~isempty(tr)
        plot(tr.x, tr.y, traceConfig{t,2}, 'LineWidth', 2)
        legEntries{end+1} = strrep(d, 'AoA', 'AoA = '); %#ok<SAGROW>
    end
end
xlabel('Flow Time'); ylabel('Lift Coefficient (C_l)')
title('Lift Coefficient (C_l) CFD Signal Trace')
if ~isempty(legEntries), legend(legEntries, 'Location', 'best'); end
ylim([0 2.5]); grid on

subplot(2,2,3)
plot(Il_cd(:,1), Il_cd(:,2), 'w-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Ku_cd(:,1), Ku_cd(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, Cd_solver, 'c-v', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Drag Coefficient (C_d)')
title('Drag Coefficient (C_d) vs AoA')
legend('Ilio','Kurtulus','Matthews','Location','best'); grid on

subplot(2,2,4)
plot(Il_clcd(:,1), Il_clcd(:,2), 'w-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
plot(Ku_clcd(:,1), Ku_clcd(:,2), 'r-s', 'LineWidth', 2, 'MarkerSize', 8)
plot(aoaValues, ClCd_calc, 'c-v', 'LineWidth', 2, 'MarkerSize', 8)
xlabel('Angle of Attack (degrees)'); ylabel('Aerodynamic Efficiency (C_l/C_d)')
title('Aerodynamic Efficiency (C_l/C_d) vs AoA')
legend('Ilio','Kurtulus','Matthews','Location','best'); grid on

saveas(gcf,'20260624_SM_dataprocess_comparison.png','png') %Rename me!

%% Plot of Reynolds' Stresses in the Domain and Wake Region
figure('Position', [300 200 1200 900])

subplot(2,3,1)
plot(aoaValues, uu_calc, 'c-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Component U (uu) [Pa]')
title('Reynolds Stress (uu) vs AoA'); legend('Matthews','Location','best'); grid on

subplot(2,3,2)
plot(aoaValues, vv_calc, 'c-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Component V (vv) [Pa]')
title('Reynolds Stress (vv) vs AoA'); legend('Matthews','Location','northwest'); grid on

subplot(2,3,3)
plot(aoaValues, uv_calc, 'c-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Component U in V (uv) [Pa]')
title('Reynolds Stress (uv) vs AoA'); legend('Matthews','Location','northwest'); grid on

subplot(2,3,4)
plot(aoaValues, uuwake_calc, 'r-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Component U (uu) [Pa]')
title('Reynolds Stress in the Wake (uu) vs AoA'); legend('Matthews','Location','northwest'); grid on

subplot(2,3,5)
plot(aoaValues, vvwake_calc, 'r-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees) [Pa]'); ylabel('Component V (vv)')
title('Reynolds Stress in the Wake (vv) vs AoA'); legend('Matthews','Location','northwest'); grid on

subplot(2,3,6)
plot(aoaValues, uvwake_calc, 'r-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Component U in V (uv) [Pa]')
title('Reynolds Stress in the Wake (uv) vs AoA'); legend('Matthews','Location','southwest'); grid on

saveas(gcf,'20260624_SM_reynolds_stress.png','png') %Rename me!

%% Wake Average Pressure Analysis
figure('Position', [200 200 1200 900])

subplot(1,2,1)
plot(aoaValues, StaticPressure_wake, 'c-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Static Pressure (Pa)')
title('Wake Static Pressure (Pa) vs AoA'); legend('Static Pressure','Location','northeast'); grid on

subplot(1,2,2)
plot(aoaValues, TotalPressure_wake, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Total Pressure (Pa)')
title('Wake Total Pressure (Pa) vs AoA'); legend('Total Pressure','Location','northeast'); grid on

saveas(gcf,'20260624_SM_wake_pressure.png','png') %Rename me!

%% Wake Average Velocity Analysis
figure('Position', [400 200 1500 900])

subplot(1,3,1)
plot(aoaValues, VelocityX_wake, 'c-v', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Wake X Velocity (m/s)')
title('Wake X Velocity vs AoA'); legend('X Velocity','Location','northeast'); grid on

subplot(1,3,2)
plot(aoaValues, VelocityY_wake, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Wake Y Velocity (m/s)')
title('Wake Y Velocity vs AoA'); legend('Y Velocity','Location','northwest'); grid on

subplot(1,3,3)
plot(aoaValues, VelocityMag_wake, 'g-s', 'LineWidth', 2, 'MarkerSize', 8); hold on
xlabel('Angle of Attack (degrees)'); ylabel('Wake Velocity Magnitude (m/s)')
title('Wake Velocity Magnitude vs AoA'); legend('Velocity Magnitude','Location','northeast'); grid on

saveas(gcf,'20260624_SM_wake_velocity.png','png') %Rename me!

%% Calculate overall difference between Matthews and Ilio data
Il_cl_interp   = interp1(Il_cl(:,1), Il_cl(:,2), aoaValues, 'linear', 'extrap');
Il_cd_interp   = interp1(Il_cd(:,1), Il_cd(:,2), aoaValues, 'linear', 'extrap');
Il_clcd_interp = interp1(Il_clcd(:,1), Il_clcd(:,2), aoaValues, 'linear', 'extrap');

Cl_percentDiff = ((Cl_solver' - Il_cl_interp) ./ Il_cl_interp) * 100;
Cl_percentDiff(1) = NaN;    % Omit first value close to zero
Cd_percentDiff = ((Cd_solver' - Il_cd_interp) ./ Il_cd_interp) * 100;
ClCd_percentDiff = ((ClCd_calc' - Il_clcd_interp) ./ Il_clcd_interp) * 100;
ClCd_percentDiff(1) = NaN;  % Omit first value divided by NaN

Cl_MAPD   = mean(abs(Cl_percentDiff), 'omitnan');
Cd_MAPD   = mean(abs(Cd_percentDiff), 'omitnan');
ClCd_MAPD = mean(abs(ClCd_percentDiff), 'omitnan');

Cl_RMSPD   = sqrt(mean(Cl_percentDiff.^2, 'omitnan'));
Cd_RMSPD   = sqrt(mean(Cd_percentDiff.^2, 'omitnan'));
ClCd_RMSPD = sqrt(mean(ClCd_percentDiff.^2, 'omitnan'));

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


%% ======================================================================= %%
%%                            LOCAL FUNCTIONS                              %%
%% ======================================================================= %%

function map = buildFileMap(dirPath)
% Build a containers.Map from normalized name -> full file path for every
% .out file in dirPath. Normalization strips case, separators, and 'rfile'.
    map = containers.Map('KeyType','char','ValueType','char');
    listing = dir(fullfile(dirPath, '*.out'));
    for k = 1:numel(listing)
        [~, stem] = fileparts(listing(k).name);
        key = normName(stem, true);   % strip the 'rfile' token for filenames
        if isempty(key), continue; end
        if isKey(map, key)
            warning('Duplicate match for key "%s" in %s (keeping first: %s).', ...
                key, dirPath, map(key));
        else
            map(key) = fullfile(dirPath, listing(k).name);
        end
    end
end

function fp = findReport(map, candidates)
% Return the path of the first candidate key present in the map, else ''.
    fp = '';
    for c = 1:numel(candidates)
        key = normName(candidates{c}, true);
        if isKey(map, key)
            fp = map(key);
            return;
        end
    end
end

function R = readReport(fp)
% Read a Fluent .out report robustly:
%   - auto-detect the number of header lines (first all-numeric line),
%   - read the column NAMES from the header and pick the value column by name,
%   - return R.y (value) and R.x (abscissa: flow-time if present, else index).
    raw = fileread(fp);
    lines = regexp(raw, '\r\n|\r|\n', 'split');

    % Find first line that is purely numeric -> data starts there.
    firstData = 0;
    for k = 1:numel(lines)
        ln = strtrim(lines{k});
        if isempty(ln), continue; end
        toks = regexp(ln, '\s+', 'split');
        vals = str2double(toks);
        if ~isempty(vals) && all(~isnan(vals))
            firstData = k;
            break;
        end
    end
    if firstData == 0
        error('No numeric data found in %s', fp);
    end
    nHead = firstData - 1;

    M = readmatrix(fp, 'FileType', 'text', 'NumHeaderLines', nHead);
    ncols = size(M, 2);

    % Parse column names: choose the header line whose count of quoted tokens
    % equals the number of data columns (skips the title line).
    names = {};
    for k = 1:nHead
        tk = regexp(lines{k}, '"([^"]*)"', 'tokens');
        if numel(tk) == ncols
            names = cellfun(@(x) x{1}, tk, 'UniformOutput', false);
        end
    end

    % Decide which column is the value, and which is the abscissa.
    valCol = []; xTimeCol = []; xIterCol = [];
    if ~isempty(names)
        for c = 1:numel(names)
            n = normName(names{c}, false);
            if contains(n, 'flowtime') || strcmp(n, 'time')
                xTimeCol = c;
            elseif contains(n, 'timestep') || contains(n, 'iteration') ...
                    || strcmp(n, 'iter') || strcmp(n, 'step')
                xIterCol = c;
            elseif isempty(valCol)
                valCol = c;   % first non-abscissa column is the value
            end
        end
    end

    % Fallbacks if names were unavailable or ambiguous.
    if isempty(valCol)
        if ncols >= 2, valCol = 2; else, valCol = 1; end
    end
    if isempty(xTimeCol) && isempty(xIterCol)
        if ncols >= 3, xTimeCol = 3; else, xIterCol = 1; end
    end

    R.y = M(:, valCol);
    if ~isempty(xTimeCol)
        R.x = M(:, xTimeCol);
    else
        R.x = M(:, xIterCol);
    end
    R.names = names;
end

function key = normName(name, stripReportToken)
% Normalize a name: lowercase, drop the extension, remove all non-alphanumeric
% characters. If stripReportToken is true, also remove the 'rfile' token so
% 'drag-rfile' and 'drag' collapse to the same key.
    if nargin < 2, stripReportToken = false; end
    [~, stem] = fileparts(name);      % harmless if there is no extension
    if isempty(stem), stem = name; end
    key = lower(stem);
    key = regexprep(key, '[^a-z0-9]', '');   % strip separators/punctuation
    if stripReportToken
        key = erase(key, 'rfile');
    end
end

function a = tailAvg(y, n)
% Mean of the trailing n samples (or all of them if fewer than n exist).
    y = y(:);
    if numel(y) > n
        y = y(end-n+1:end);
    end
    a = mean(y, 'omitnan');
end