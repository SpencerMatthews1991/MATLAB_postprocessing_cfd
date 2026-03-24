%% FLUENT JOURNAL FILE GENERATOR FOR ANGLE OF ATTACK SWEEPS
% This script reads AoA lookup data from CSV and generates Fluent journal
% files for each AoA case folder with updated boundary conditions.

% ---------------- HARDCODED CONFIGURATION ---------------- %
ROOT_DIR = pwd;  % Uses current working directory
CSV_FILE = '20251126_SM_ForceTransform_Ver1.csv';
ITERATIONS = 100;  % Number of iterations for Fluent solve

% ---------------- REFERENCE VALUES CONFIGURATION ---------------- %
% Values must match your Fluent Reference Values panel (Setup > Reference Values)
REF.compute_from  = 'inlet_curves';  % Zone to compute from
REF.area          = 1;               % Area [m^2]
REF.density       = 1;               % Density [kg/m^3]
REF.depth         = 1;               % Depth [m]
REF.enthalpy      = 0;               % Enthalpy [J/kg]
REF.length        = 1;               % Length [m]
REF.pressure      = 101325;          % Pressure [Pa]
REF.temperature   = 288.16;          % Temperature [K]
REF.velocity      = 1;               % Velocity [m/s]
REF.viscosity     = 0.001;           % Viscosity [kg/(m s)]
REF.gamma         = 1.4;             % Ratio of Specific Heats
REF.yplus         = 300;             % Yplus for Heat Transfer Coefficient
REF.zone          = 'fluid';         % Reference Zone

%% MAIN SCRIPT
fprintf('Processing directory: %s\n', ROOT_DIR);
fprintf('Using CSV file: %s\n', CSV_FILE);
fprintf('Iterations: %d\n', ITERATIONS);
fprintf('%s\n', repmat('-', 1, 50));

% Debug: Check if CSV file exists
csv_full_path = fullfile(ROOT_DIR, CSV_FILE);
fprintf('Looking for CSV at: %s\n', csv_full_path);
fprintf('CSV file exists: %d\n', exist(csv_full_path, 'file'));

% List files in directory for debugging
files = dir(ROOT_DIR);
files = files(~[files.isdir]);  % Remove directories
fprintf('Files in directory (first 5): ');
for i = 1:min(5, length(files))
    fprintf('%s ', files(i).name);
end
fprintf('\n%s\n', repmat('-', 1, 50));

% Read lookup table from CSV
lookup = read_lookup(csv_full_path);

% Get all directories in ROOT_DIR
folders = dir(ROOT_DIR);
folders = folders([folders.isdir]);  % Keep only directories
folder_names = {folders.name};

% Filter for AoA folders and sort them
aoa_folders = folder_names(startsWith(folder_names, 'AoA'));
aoa_folders = sort(aoa_folders);  % Sort alphabetically/numerically

% Process each AoA folder
for i = 1:length(aoa_folders)
    folder = aoa_folders{i};
    
    % Extract AoA value from folder name
    aoa_str = strrep(folder, 'AoA', '');
    aoa = str2double(aoa_str);
    
    if isnan(aoa)
        fprintf('Skipping folder with invalid AoA format: %s\n', folder);
        continue;
    end
    
    % Check if lookup data exists for this AoA
    if ~isKey(lookup, aoa)
        fprintf('No lookup data for AoA %.1f\n', aoa);
        continue;
    end
    
    data = lookup(aoa);
    folder_path = fullfile(ROOT_DIR, folder);
    
    % Find .cas.h5 file in folder
    case_file = find_case_file(folder_path);
    if isempty(case_file)
        fprintf('No .cas.h5 file found in %s\n', folder);
        continue;
    end
    
    % Generate journal file
    generate_journal(folder_path, aoa, case_file, data, ITERATIONS, REF);
end

fprintf('%s\n', repmat('-', 1, 50));
fprintf('Done!\n');

%% FUNCTION: Read lookup table from CSV
function lookup = read_lookup(csv_file)
    % Read CSV file
    data = readtable(csv_file);
    
    % Create containers.Map (dictionary) for lookup
    lookup = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    % Process each row
    for i = 1:height(data)
        aoa = data.AoA(i);
        
        % Create structure with all necessary data
        entry = struct();
        entry.Ux = num2str(data.Ux(i));  % Convert to string
        entry.Uy = num2str(data.Uy(i));
        entry.Uz = '0';
        entry.drag = {num2str(data.drag_x(i)), num2str(data.drag_y(i)), '0'};
        entry.lift = {num2str(data.lift_x(i)), num2str(data.lift_y(i)), '0'};
        
        lookup(aoa) = entry;
    end
end

%% FUNCTION: Find .cas.h5 file in folder
function case_file = find_case_file(folder_path)
    files = dir(fullfile(folder_path, '*.cas.h5'));
    
    if isempty(files)
        case_file = '';
    else
        case_file = files(1).name;  % Return first match
    end
end

%% FUNCTION: Generate Fluent journal file
function generate_journal(folder_path, aoa, case_file, data, iterations, ref)
    journal_path = fullfile(folder_path, sprintf('update_AoA%d.jou', round(aoa)));
    
    % Use relative path - just the filename since journal is in same folder
    case_file_path = case_file;
    
    % Generate new date-stamped filenames
    today = datestr(now, 'yyyymmdd');
    aoa_int = round(aoa);
    base_name = sprintf('%s_SM_mediumNACA0012_NonPorous_AoA%d_Ver5', today, aoa_int);
    
    case_out = sprintf('%s.cas.h5', base_name);
    data_out = sprintf('%s.dat.h5', base_name);
    export_name = sprintf('./Exports_AoA%d/%s_SM_mediumNACA0012_NonPorous_AoA%d_Ver5', ...
                          aoa_int, today, aoa_int);
    autosave_name = sprintf('./%s_SM_mediumNACA0012_NonPorous_AoA%d_Ver5', today, aoa_int);
    
    % Open file for writing
    fid = fopen(journal_path, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', journal_path);
    end
    
    % ------------------------------------------------------------------ %
    % Read the case file
    % ------------------------------------------------------------------ %
    fprintf(fid, '/file/read-case %s\n', case_file_path);
    fprintf(fid, '\n');

    % ------------------------------------------------------------------ %
    % SET REFERENCE VALUES (TUI)
    % Equivalent to: Setup > Reference Values panel in the GUI
    %
    % TUI path: /report/reference-values/
    %   compute          – populates from a boundary zone
    %   area             – reference area  [m^2]
    %   density          – reference density [kg/m^3]
    %   depth            – reference depth [m]
    %   enthalpy         – reference enthalpy [J/kg]
    %   length           – reference length [m]
    %   pressure         – reference pressure [Pa]
    %   temperature      – reference temperature [K]
    %   velocity         – reference velocity [m/s]
    %   viscosity        – reference viscosity [kg/(m s)]
    %   ratio-of-sp-heats – ratio of specific heats (gamma)
    %   yplus-for-htc    – y+ for heat transfer coefficient
    %
    % Note: 'compute' sets most values automatically from the chosen zone.
    % The individual commands below then override each field explicitly,
    % ensuring the values match the GUI screenshot regardless of what
    % Fluent would compute automatically.
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Reference Values ---\n');

    % Compute from inlet zone first (populates defaults)
    fprintf(fid, '/report/reference-values/compute %s\n', ref.compute_from);

    % Override each field explicitly with configured values
    fprintf(fid, '/report/reference-values/area %g\n',        ref.area);
    fprintf(fid, '/report/reference-values/density %g\n',     ref.density);
    fprintf(fid, '/report/reference-values/depth %g\n',       ref.depth);
    fprintf(fid, '/report/reference-values/enthalpy %g\n',    ref.enthalpy);
    fprintf(fid, '/report/reference-values/length %g\n',      ref.length);
    fprintf(fid, '/report/reference-values/pressure %g\n',    ref.pressure);
    fprintf(fid, '/report/reference-values/temperature %g\n', ref.temperature);
    fprintf(fid, '/report/reference-values/velocity %g\n',    ref.velocity);
    fprintf(fid, '/report/reference-values/viscosity %g\n',   ref.viscosity);
    fprintf(fid, '/report/reference-values/ratio-of-sp-heats %g\n', ref.gamma);
    fprintf(fid, '/report/reference-values/yplus-for-htc %g\n',     ref.yplus);
    fprintf(fid, '\n');

    % ------------------------------------------------------------------ %
    % Set velocity inlet boundary condition (2D case)
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Velocity Inlet Boundary Condition ---\n');
    fprintf(fid, '/define/boundary-conditions/velocity-inlet\n');
    fprintf(fid, 'inlet_curves\n');    % Zone name
    fprintf(fid, 'no\n');              % Magnitude and Direction? no
    fprintf(fid, 'yes\n');             % Components? yes
    fprintf(fid, 'yes\n');             % Absolute reference frame? yes
    fprintf(fid, 'no\n');              % Profile for Supersonic/Initial Gauge Pressure? no
    fprintf(fid, '0\n');               % Supersonic/Initial Gauge Pressure value
    fprintf(fid, 'no\n');              % Profile for X-velocity? no
    fprintf(fid, '%s\n', data.Ux);     % X-velocity value
    fprintf(fid, 'no\n');              % Profile for Y-velocity? no
    fprintf(fid, '%s\n', data.Uy);     % Y-velocity value
    fprintf(fid, '\n');                % Exit velocity-inlet menu
    
    % ------------------------------------------------------------------ %
    % Set force vectors in report definitions (2D – X and Y only)
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Report Definition Force Vectors ---\n');
    fprintf(fid, '/solve/report-definitions/edit lift force-vector %s %s q\n', ...
            data.lift{1}, data.lift{2});
    fprintf(fid, '/solve/report-definitions/edit drag force-vector %s %s q\n', ...
            data.drag{1}, data.drag{2});
    fprintf(fid, '/solve/report-definitions/edit lift_coefficient force-vector %s %s q\n', ...
            data.lift{1}, data.lift{2});
    fprintf(fid, '/solve/report-definitions/edit drag_coefficient force-vector %s %s q\n', ...
            data.drag{1}, data.drag{2});
    fprintf(fid, '\n');

    % ------------------------------------------------------------------ %
    % Initialize and solve
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Initialize and Solve ---\n');
    fprintf(fid, '/solve/initialize/initialize-flow\n');
    
    % ------------------------------------------------------------------ %
    % Write updated case and data files
    % ------------------------------------------------------------------ %
    fprintf(fid, '/file/write-case-data %s yes\n', case_out);
    
    % ------------------------------------------------------------------ %
    % Exit Fluent
    % ------------------------------------------------------------------ %
    fprintf(fid, '/exit yes\n');
    
    % Close file
    fclose(fid);
    
    % Print confirmation
    fprintf('Journal file created: %s\n', journal_path);
    fprintf('  Output files: %s, %s\n', case_out, data_out);
    fprintf('  Export prefix: %s\n', export_name);
    fprintf('  Autosave prefix: %s\n', autosave_name);
end
