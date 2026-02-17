%% Fluent Live Residual Monitor
% Reads an ANSYS Fluent transcript file and plots continuity,
% x-velocity, and y-velocity residuals in real-time as the
% simulation progresses.
%
% Usage: Run this script while your Fluent simulation is running.
%        Update 'transcriptFile' to point to your transcript path.

clear; clc; close all;

%% --- USER SETTINGS ---------------------------------------------------
searchDir = 'C:\Users\Spencer.Matthews.397\Desktop\Porous _Airfoil _Reseach\20260206_SM_test_case_AoA16\AoA16\';
filePattern = '*.trn';

% Find file matching wildcard
matches = dir(fullfile(searchDir, filePattern));
if isempty(matches)
    error('No file matching "%s" found in:\n%s', filePattern, searchDir);
elseif numel(matches) > 1
    [~, idx] = max([matches.datenum]);
    matches = matches(idx);
    fprintf('Multiple matches found – using most recent: %s\n', matches.name);
end
transcriptFile = fullfile(matches.folder, matches.name);

refreshInterval = 2;        % seconds between file re-reads
maxPoints       = 0;     % max data points to keep on plot (0 = unlimited)
% ----------------------------------------------------------------------

%% Set up figure
fig = figure('Name', 'Fluent Live Residuals', ...
             'NumberTitle', 'off', ...
             'Color', 'white', ...
             'Position', [100 100 900 550]);

ax = axes('Parent', fig);
semilogy(ax, NaN, NaN);     % initialise log-scale axes
hold(ax, 'on');
grid(ax, 'on');
grid(ax, 'minor');

% Line handles – created empty so we can update them in the loop
hCont  = semilogy(ax, NaN, NaN, '-o', 'Color', [0.85 0.15 0.10], ...
                  'MarkerSize', 3, 'DisplayName', 'Continuity');
hVelX  = semilogy(ax, NaN, NaN, '-s', 'Color', [0.10 0.45 0.85], ...
                  'MarkerSize', 3, 'DisplayName', 'x-velocity');
hVelY  = semilogy(ax, NaN, NaN, '-^', 'Color', [0.10 0.70 0.30], ...
                  'MarkerSize', 3, 'DisplayName', 'y-velocity');

xlabel(ax, 'Time Step',   'FontSize', 12);
ylabel(ax, 'Residual',    'FontSize', 12);
title(ax,  'ANSYS Fluent Residuals (live)', 'FontSize', 14, 'FontWeight', 'bold');
legend(ax, 'Location', 'northeast', 'FontSize', 11);
ax.YScale = 'log';

% Annotation for current flow-time
hTime = annotation(fig, 'textbox', [0.13 0.01 0.5 0.05], ...
                   'String', 'Flow time: --', ...
                   'EdgeColor', 'none', ...
                   'FontSize', 10, 'Color', [0.4 0.4 0.4]);

% Status indicator in title bar
fprintf('Monitoring: %s\n', transcriptFile);
fprintf('Press Ctrl+C in the Command Window to stop.\n\n');

%% Regex pattern
% Matches lines like:
%    100  1.6286e-04  3.0675e-08  2.8995e-08  0:00:00    0
% (the INNER sub-iteration line, NOT the outer "step" line)
pattern = ['^\s*(\d+)\s+'                    ... % iter number
           '([0-9]+\.[0-9]+e[+-][0-9]+)\s+'  ... % continuity
           '([0-9]+\.[0-9]+e[+-][0-9]+)\s+'  ... % x-velocity
           '([0-9]+\.[0-9]+e[+-][0-9]+)'];        % y-velocity

stepPattern  = 'step\s+flow-time\s+delta-time';
ftimePattern = '^\s*(\d+)\s+([\d.e+-]+)\s+([\d.e+-]+)';

%% Live loop
prevLines = 0;      % track how many lines we have already processed
steps  = [];
cont   = [];
velX   = [];
velY   = [];
ftimes = {};

while ishandle(fig)
    % ---- Read file ---------------------------------------------------
    fid = fopen(transcriptFile, 'r');
    if fid == -1
        title(ax, sprintf('Waiting for file: %s', transcriptFile), ...
              'FontSize', 12, 'Color', 'r');
        pause(refreshInterval);
        continue
    end
    raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = raw{1};

    % Only process new lines since last read
    newLines = lines(prevLines+1 : end);
    prevLines = numel(lines);

    % ---- Parse new lines ---------------------------------------------
    awaitingStep = false;
    for k = 1 : numel(newLines)
        ln = strtrim(newLines{k});

        % Detect the "step  flow-time  delta-time" header
        if ~isempty(regexp(ln, stepPattern, 'once'))
            awaitingStep = true;
            continue
        end

        % The line immediately after the header gives step number & flow-time
        if awaitingStep
            tok = regexp(ln, ftimePattern, 'tokens', 'once');
            if ~isempty(tok)
                steps(end+1) = str2double(tok{1}); %#ok<AGROW>
                ftimes{end+1} = tok{2};            %#ok<AGROW>
            end
            awaitingStep = false;
            continue
        end

        % Residual line (inner iteration with 4 scientific numbers)
        tok = regexp(ln, pattern, 'tokens', 'once');
        if numel(tok) == 4
            cont(end+1) = str2double(tok{2}); %#ok<AGROW>
            velX(end+1) = str2double(tok{3}); %#ok<AGROW>
            velY(end+1) = str2double(tok{4}); %#ok<AGROW>
        end
    end

    % ---- Update plot -------------------------------------------------
    n = min([numel(cont), numel(velX), numel(velY)]);
    if n > 0
        idx = max(1, n - maxPoints + 1) : n;   % window if maxPoints set
        if maxPoints == 0, idx = 1:n; end

        set(hCont, 'XData', idx, 'YData', cont(idx));
        set(hVelX, 'XData', idx, 'YData', velX(idx));
        set(hVelY, 'XData', idx, 'YData', velY(idx));

        % Update flow-time label
        if ~isempty(ftimes)
            set(hTime, 'String', sprintf('Flow time: %s s   |   Step: %d', ...
                ftimes{end}, numel(steps)));
        end

        % Auto-scale y limits with a small margin
        allVals = [cont(idx), velX(idx), velY(idx)];
        allVals = allVals(allVals > 0);
        if ~isempty(allVals)
            ylo = 10^(floor(log10(min(allVals))) - 0.5);
            yhi = 10^(ceil( log10(max(allVals))) + 0.5);
            ylim(ax, [ylo, yhi]);
        end
        xlim(ax, [max(0, idx(1)-1), idx(end)+1]);

        drawnow limitrate
    end

    pause(refreshInterval);
end

fprintf('Figure closed – monitoring stopped.\n')