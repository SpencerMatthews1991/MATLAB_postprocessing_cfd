%% FLUENT JOURNAL FILE GENERATOR FOR ANGLE OF ATTACK SWEEPS
% This script reads AoA lookup data from CSV and generates Fluent journal
% files for each AoA case folder with updated boundary conditions.
% Ver4: Added wake rake generation - multiple parallel rakes perpendicular
%       to the wake centreline, adjusted per AoA.
% Ver6: Fixed rake Type dropdown command for Fluent 2025 R1 (cx-set-list-selections).

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

% ---------------- WAKE RAKE CONFIGURATION ---------------- %
% All distances in chord lengths (chord = 1.0 m assumed)
RAKE.x_start      = 1.1;    % Downstream start of rake from origin [c]
RAKE.x_end        = 30.1;    % Downstream end of rake [c]
RAKE.n_points     = 120;     % Number of points per rake
RAKE.n_rakes      = 3;      % Number of parallel rakes (odd number recommended)
RAKE.delta        = 0.15;    % Perpendicular spacing between rakes [c]
RAKE.fluid_zone   = 'fluid'; % Fluid zone name in Fluent

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
    generate_journal(folder_path, aoa, case_file, data, ITERATIONS, REF, RAKE);
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

%% FUNCTION: Compute rake endpoints for a given AoA
function rakes = compute_rake_coords(aoa, rake_cfg)
    % Computes start and end coordinates for each rake in the array,
    % rotated to align with the wake direction at the given AoA.
    %
    % Wake direction is parallel to the freestream (along AoA vector).
    % Rakes are oriented perpendicular to the wake, offset by delta.
    %
    % Inputs:
    %   aoa       - angle of attack in degrees
    %   rake_cfg  - struct with fields: x_start, x_end, n_points,
    %               n_rakes, delta
    %
    % Output:
    %   rakes     - struct array with fields: x1, y1, x2, y2, name

    aoa_rad = deg2rad(aoa);

    % Unit vector along wake (freestream direction)
    along = [cos(aoa_rad);  sin(aoa_rad)];

    % Unit vector perpendicular to wake (normal direction in 2D)
    perp  = [-sin(aoa_rad); cos(aoa_rad)];

    % Start and end points along wake centreline
    p1_centre = rake_cfg.x_start * along;
    p2_centre = rake_cfg.x_end   * along;

    % Offsets: symmetric about centreline, e.g. n_rakes=3 -> [-1, 0, 1]*delta
    half = floor(rake_cfg.n_rakes / 2);
    offsets = (-half : half) * rake_cfg.delta;  % length = n_rakes

    rakes(rake_cfg.n_rakes) = struct('x1',[],'y1',[],'x2',[],'y2',[],'name','');

    for k = 1:rake_cfg.n_rakes
        offset = offsets(k);
        p1 = p1_centre + offset * perp;
        p2 = p2_centre + offset * perp;

        rakes(k).x1   = p1(1);
        rakes(k).y1   = p1(2);
        rakes(k).x2   = p2(1);
        rakes(k).y2   = p2(2);
        rakes(k).name = sprintf('wake_rake_%d', k);
    end
end

%% FUNCTION: Generate Fluent journal file
function generate_journal(folder_path, aoa, case_file, data, iterations, ref, rake_cfg)
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
    
    % Compute rake coordinates for this AoA
    rakes = compute_rake_coords(aoa, rake_cfg);

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
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Reference Values ---\n');
    % fprintf(fid, '/report/reference-values/compute %s\n', ref.compute_from);  % Removed - causes error in 2025 R1
    fprintf(fid, '/report/reference-values/area %g\n',        ref.area);
    fprintf(fid, '/report/reference-values/density %g\n',     ref.density);
    fprintf(fid, '/report/reference-values/depth %g\n',       ref.depth);
    fprintf(fid, '/report/reference-values/enthalpy %g\n',    ref.enthalpy);
    fprintf(fid, '/report/reference-values/length %g\n',      ref.length);
    fprintf(fid, '/report/reference-values/pressure %g\n',    ref.pressure);
    fprintf(fid, '/report/reference-values/temperature %g\n', ref.temperature);
    fprintf(fid, '/report/reference-values/velocity %g\n',    ref.velocity);
    fprintf(fid, '/report/reference-values/viscosity %g\n',   ref.viscosity);
    % fprintf(fid, '/report/reference-values/ratio-of-sp-heats %g\n', ref.gamma);  % Removed - invalid in 2025 R1 laminar case
    % fprintf(fid, '/report/reference-values/yplus-for-htc %g\n',     ref.yplus);  % Removed - not relevant for laminar incompressible
    fprintf(fid, '\n');

    % ------------------------------------------------------------------ %
    % Set velocity inlet boundary condition (2D case)
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Velocity Inlet Boundary Condition ---\n');
    fprintf(fid, '/define/boundary-conditions/velocity-inlet\n');
    fprintf(fid, 'inlet_curves\n');
    fprintf(fid, 'no\n');
    fprintf(fid, 'yes\n');
    fprintf(fid, 'yes\n');
    fprintf(fid, 'no\n');
    fprintf(fid, '0\n');
    fprintf(fid, 'no\n');
    fprintf(fid, '%s\n', data.Ux);
    fprintf(fid, 'no\n');
    fprintf(fid, '%s\n', data.Uy);
    fprintf(fid, '\n');
    
    % ------------------------------------------------------------------ %
    % Set force vectors in report definitions (2D - X and Y only)
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
    % CREATE wake rakes adjusted for this AoA
    %
    % Uses GUI commands (cx-gui-do) as recorded in Fluent 2025 R1.
    % The Line/Rake Surface dialog expects:
    %   - Surface name
    %   - x0 (start x), y0 (start y)
    %   - x1 (end x),   y1 (end y)
    %   - Number of points
    %
    % Note: RealEntry3 and RealEntry6 are z-coordinates (default 0 in 2D).
    % Rakes are parallel lines offset perpendicular to the wake centreline.
    % ------------------------------------------------------------------ %
    fprintf(fid, '; --- Wake Rake Surfaces (AoA = %g deg) ---\n', aoa);
    fprintf(fid, '; along=[cos(aoa),sin(aoa)], perp=[-sin(aoa),cos(aoa)]\n');
    fprintf(fid, '; x_start=%.3fc, x_end=%.3fc, %d pts/rake, %d rakes, delta=%.3fc\n', ...
            rake_cfg.x_start, rake_cfg.x_end, rake_cfg.n_points, ...
            rake_cfg.n_rakes, rake_cfg.delta);
    fprintf(fid, '\n');

    for k = 1:rake_cfg.n_rakes
        % Comment line with coords for human readability
        fprintf(fid, '; Rake %d  x0=%.4f y0=%.4f  x1=%.4f y1=%.4f\n', ...
                k, rakes(k).x1, rakes(k).y1, rakes(k).x2, rakes(k).y2);

        % Open the Line/Rake dialog via the Results ribbon
        fprintf(fid, '(cx-gui-do cx-activate-item "Ribbon*Frame1*Frame6(Results)*Table1*Table3(Surface)*PushButton1(Create)")\n');
        fprintf(fid, '(cx-gui-do cx-activate-item "MenuBar*PopupMenuCreate*Line/Rake...")\n');

        % Switch Type dropdown from Line to Rake
        fprintf(fid, '(cx-gui-do cx-set-list-selections "Line/Rake Surface*Table1*Table2(Options)*DropDownList2(Type)" ''("Rake"))\n');

        % Set surface name
        fprintf(fid, ['(cx-gui-do cx-set-text-entry ' ...
            '"Line/Rake Surface*Table1*TextEntry1(New Surface Name)" "%s")\n'], ...
            rakes(k).name);

        % Set start point x0, y0
        fprintf(fid, ['(cx-gui-do cx-set-real-entry-list ' ...
            '"Line/Rake Surface*Table1*Frame5(End Points)*Table1*RealEntry1(x0)" ' ...
            '''( %g))\n'], rakes(k).x1);
        fprintf(fid, ['(cx-gui-do cx-set-real-entry-list ' ...
            '"Line/Rake Surface*Table1*Frame5(End Points)*Table1*RealEntry2(y0)" ' ...
            '''( %g))\n'], rakes(k).y1);

        % Set end point x1, y1
        fprintf(fid, ['(cx-gui-do cx-set-real-entry-list ' ...
            '"Line/Rake Surface*Table1*Frame5(End Points)*Table1*RealEntry4(x1)" ' ...
            '''( %g))\n'], rakes(k).x2);
        fprintf(fid, ['(cx-gui-do cx-set-real-entry-list ' ...
            '"Line/Rake Surface*Table1*Frame5(End Points)*Table1*RealEntry5(y1)" ' ...
            '''( %g))\n'], rakes(k).y2);

        % Set number of points
        fprintf(fid, ['(cx-gui-do cx-set-integer-entry ' ...
            '"Line/Rake Surface*Table1*IntegerEntry1(Number of Points)" %d)\n'], ...
            rake_cfg.n_points);

        % Confirm dialog
        fprintf(fid, '(cx-gui-do cx-activate-item "Line/Rake Surface*PanelButtons*PushButton1(OK)")\n');
        fprintf(fid, '\n');
    end


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
    fprintf('  AoA = %g deg | %d rakes x %d pts | delta = %.3fc\n', ...
            aoa, rake_cfg.n_rakes, rake_cfg.n_points, rake_cfg.delta);
    fprintf('  Output files: %s, %s\n', case_out, data_out);
    fprintf('  Export prefix: %s\n', export_name);
    fprintf('  Autosave prefix: %s\n', autosave_name);
end