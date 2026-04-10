%% GCI script, version 1
% Initialization
clear all
close all
clc

%% Find directory
baseDir = fullfile(fileparts(mfilename('fullpath')), 'Code_Testing');

%% Mesh information from user inputs (spatial)
Mesh_supercoarse = 12288;
Mesh_mediumcoarse = 25284;
Mesh_coarse = 51527;
Mesh_medium = 206108;
Mesh_fine   = 824432;

%% Mesh information (temporal)
dt_supercoarse = 0.001;      %sec
dt_mediumcoarse = 0.0006667; %sec
dt_coarse = 0.0005;          %sec
dt_medium = 0.00025;         %sec
dt_fine = 0.0000125;         %sec

%% Import GCI mesh data from files
supercoarse = readtable(fullfile(baseDir, '20260324_AoA_16_GCI_supercoarse.csv'));
mediumcoarse = readtable(fullfile(baseDir, '20260324_AoA_16_GCI_mediumcoarse.csv'));
coarse = readtable(fullfile(baseDir, '20260213_AoA_16_GCI_coarse.csv'));
medium = readtable(fullfile(baseDir, '20260213_AoA_16_GCI_medium.csv'));
fine   = readtable(fullfile(baseDir, '20260213_AoA_16_GCI_fine.csv'));

% Extract Cl, Cd, Cl/Cd (columns 5, 6, 7)
Grid1_cl   = supercoarse{1, 5};
Grid1_cd   = supercoarse{1, 6};
Grid1_clcd = supercoarse{1, 7};

Grid2_cl   = mediumcoarse{1, 5};
Grid2_cd   = mediumcoarse{1, 6};
Grid2_clcd = mediumcoarse{1, 7};

Grid3_cl   = coarse{1, 5};
Grid3_cd   = coarse{1, 6};
Grid3_clcd = coarse{1, 7};

Grid4_cl   = medium{1, 5};
Grid4_cd   = medium{1, 6};
Grid4_clcd = medium{1, 7};

Grid5_cl   = fine{1, 5};
Grid5_cd   = fine{1, 6};
Grid5_clcd = fine{1, 7};

%% Mesh spacing (2D) temporal 
h_supercoarse  = sqrt(1/Mesh_supercoarse);
h_mediumcoarse = sqrt(1/Mesh_mediumcoarse);
h_coarse       = sqrt(1/Mesh_coarse);
h_medium       = sqrt(1/Mesh_medium);
h_fine         = sqrt(1/Mesh_fine);

% Normalise all by finest mesh spacing
h_norm_supercoarse  = h_supercoarse  / h_fine;
h_norm_mediumcoarse = h_mediumcoarse / h_fine;
h_norm_coarse       = h_coarse       / h_fine;
h_norm_medium       = h_medium       / h_fine;
h_norm_fine         = h_fine         / h_fine; % = 1.0

h_array = [h_fine, h_medium, h_coarse, h_mediumcoarse, h_supercoarse];
h_norm = [h_norm_fine, h_norm_medium, h_norm_coarse, h_norm_mediumcoarse, h_norm_supercoarse];

%% Mesh spacing temporal
t_supercoarse  = sqrt(1/dt_supercoarse);
t_mediumcoarse = sqrt(1/dt_mediumcoarse);
t_coarse       = sqrt(1/dt_coarse);
t_medium       = sqrt(1/dt_medium);
t_fine         = sqrt(1/dt_fine);

% Normalise all by finest mesh spacing
t_norm_supercoarse  = t_supercoarse  / t_fine;
t_norm_mediumcoarse = t_mediumcoarse / t_fine;
t_norm_coarse       = t_coarse       / t_fine;
t_norm_medium       = t_medium       / t_fine;
t_norm_fine         = t_fine         / t_fine; % = 1.0

t_array = [t_fine, t_medium, t_coarse, t_mediumcoarse, t_supercoarse];
t_norm = [t_norm_fine, t_norm_medium, t_norm_coarse, t_norm_mediumcoarse, t_norm_supercoarse];

%% Refinement ratios
r_scmc = sqrt(Mesh_mediumcoarse / Mesh_supercoarse); 
r_mcc  = sqrt(Mesh_coarse / Mesh_mediumcoarse);       
r_cm   = sqrt(Mesh_medium / Mesh_coarse);             
r_mf   = sqrt(Mesh_fine   / Mesh_medium);

fprintf('=== Refinement Ratios ===\n');
fprintf('Supercoarse to Mediumcoarse: %.4f\n', r_scmc);
fprintf('Mediumcoarse to Coarse: %.4f\n', r_mcc);
fprintf('Coarse to Medium: %.4f\n', r_cm);
fprintf('Medium to Fine:   %.4f\n', r_mf);

%% Observed order of accuracy for each variable
p_cd   = log((Grid1_cd   - Grid2_cd)  /(Grid2_cd   - Grid3_cd)) / log(r_scmc);
p_cl   = log((Grid1_cl   - Grid2_cl)  /(Grid2_cl   - Grid3_cl)) / log(r_scmc);
p_clcd = log((Grid1_clcd - Grid2_clcd)/(Grid2_clcd - Grid3_clcd)) / log(r_scmc);

fprintf('\n=== Observed Order of Accuracy ===\n');
fprintf('p (Cd):   %.4f\n', p_cd);
fprintf('p (Cl):   %.4f\n', p_cl);
fprintf('p (Cl/Cd): %.4f\n', p_clcd);

%% GCI calculations for each variable
% Using observed p for each variable
% Cd
epsilon_cd_12 = (Grid2_cd - Grid1_cd) / Grid1_cd;
epsilon_cd_23 = (Grid3_cd - Grid2_cd) / Grid2_cd;
epsilon_cd_34 = (Grid4_cd - Grid3_cd) / Grid3_cd;
epsilon_cd_45 = (Grid5_cd - Grid4_cd) / Grid4_cd;
GCI_cd_supercoarse = (3 * abs(epsilon_cd_12)) / (r_scmc^p_cd - 1);
GCI_cd_coarse = (3 * abs(epsilon_cd_23)) / (r_mcc^p_cd - 1);
GCI_cd_medium = (3 * abs(epsilon_cd_34)) / (r_cm^p_cd - 1);
GCI_cd_fine = (3 * abs(epsilon_cd_45)) / (r_mf^p_cd - 1);
f_exact_cd   = Grid3_cd + (Grid3_cd - Grid2_cd) / (r_mcc^p_cd   - 1);

% Cl
epsilon_cl_12 = (Grid2_cl - Grid1_cl) / Grid1_cl;
epsilon_cl_23 = (Grid3_cl - Grid2_cl) / Grid2_cl;
epsilon_cl_34 = (Grid4_cl - Grid3_cl) / Grid3_cl;
epsilon_cl_45 = (Grid5_cl - Grid4_cl) / Grid4_cl;
GCI_cl_supercoarse = (3 * abs(epsilon_cl_12)) / (r_scmc^p_cl - 1);
GCI_cl_coarse = (3 * abs(epsilon_cl_23)) / (r_mcc^p_cl - 1);
GCI_cl_medium = (3 * abs(epsilon_cl_34)) / (r_cm^p_cl - 1);
GCI_cl_fine = (3 * abs(epsilon_cl_45)) / (r_mf^p_cl - 1);
f_exact_cl   = Grid3_cl + (Grid3_cl - Grid2_cl) / (r_mcc^p_cl - 1);

% Cl/Cd
epsilon_clcd_12 = (Grid2_clcd - Grid1_clcd) / Grid1_clcd;
epsilon_clcd_23 = (Grid3_clcd - Grid2_clcd) / Grid2_clcd;
epsilon_clcd_34 = (Grid4_clcd - Grid3_clcd) / Grid3_clcd;
epsilon_clcd_45 = (Grid5_clcd - Grid4_clcd) / Grid4_clcd;
GCI_clcd_supercoarse = (3 * abs(epsilon_clcd_12)) / (r_scmc^p_clcd - 1);
GCI_clcd_coarse = (3 * abs(epsilon_clcd_23)) / (r_mcc^p_clcd - 1);
GCI_clcd_medium = (3 * abs(epsilon_clcd_34)) / (r_cm^p_clcd - 1);
GCI_clcd_fine = (3 * abs(epsilon_clcd_45)) / (r_mf^p_clcd - 1);
f_exact_clcd = Grid3_clcd + (Grid3_clcd - Grid2_clcd) / (r_mcc^p_clcd - 1);

%% GCI - update asymptotic ratio to use r_mf for medium pair
asymptotic_cd = GCI_cd_coarse / (r_mcc^p_cd * GCI_cd_supercoarse);
asymptotic_cl   = GCI_cl_medium   / (r_mf^p_cl   * GCI_cl_coarse);
asymptotic_clcd = GCI_clcd_medium / (r_mf^p_clcd * GCI_clcd_coarse);

%% Print results
fprintf('\n=========== GCI Results: Cd monotonic ===========\n');
fprintf('        Grid1 Cd:      %.6f\n', Grid1_cd);
fprintf('        Grid2 Cd:      %.6f\n', Grid2_cd);
fprintf('        Grid3 Cd:      %.6f\n', Grid3_cd);
fprintf('      f_exact Cd:      %.6f\n', f_exact_cd);
fprintf(' GCI_supercoarse:      %.6f (%.4f%%)\n', GCI_cd_supercoarse, GCI_cd_supercoarse*100);
fprintf('  GCI_medcoarse:       %.6f (%.4f%%)\n', GCI_cd_coarse, GCI_cd_coarse*100);
fprintf('Asymptotic ratio:      %.4f (should be ~1.0)\n', asymptotic_cd);
fprintf('\n*** Grid4 and Grid5 excluded: non-monotonic convergence ***\n');

fprintf('\n=========== GCI Results: Cd ===========\n');
fprintf('        Grid1 Cd:      %.6f\n', Grid1_cd);
fprintf('        Grid2 Cd:      %.6f\n', Grid2_cd);
fprintf('        Grid3 Cd:      %.6f\n', Grid3_cd);
fprintf('        Grid4 Cd:      %.6f\n', Grid4_cd);
fprintf('        Grid5 Cd:      %.6f\n', Grid5_cd);
fprintf('      f_exact Cd:      %.6f\n', f_exact_cd);
fprintf(' GCI_supercoarse:      %.6f (%.4f%%)\n', GCI_cd_supercoarse, GCI_cd_supercoarse*100);
fprintf('      GCI_coarse:      %.6f (%.4f%%)\n', GCI_cd_coarse, GCI_cd_coarse*100);
fprintf('      GCI_medium:      %.6f (%.4f%%)\n', GCI_cd_medium, GCI_cd_medium*100);
fprintf('        GCI_fine:      %.6f (%.4f%%)\n', GCI_cd_fine, GCI_cd_fine*100);
fprintf('Asymptotic ratio:      %.4f (should be ~1.0)\n', asymptotic_cd);

fprintf('\n========== GCI Results: Cl ===========\n');
fprintf('        Grid1 Cl:      %.6f\n', Grid1_cl);
fprintf('        Grid2 Cl:      %.6f\n', Grid2_cl);
fprintf('        Grid3 Cl:      %.6f\n', Grid3_cl);
fprintf('        Grid4 Cl:      %.6f\n', Grid4_cl);
fprintf('        Grid5 Cl:      %.6f\n', Grid5_cl);
fprintf('      f_exact Cl:      %.6f\n', f_exact_cl);
fprintf(' GCI_supercoarse:      %.6f (%.4f%%)\n', GCI_cl_supercoarse, GCI_cl_supercoarse*100);
fprintf('      GCI_coarse:      %.6f (%.4f%%)\n', GCI_cl_coarse, GCI_cl_coarse*100);
fprintf('      GCI_medium:      %.6f (%.4f%%)\n', GCI_cl_medium, GCI_cl_medium*100);
fprintf('        GCI_fine:      %.6f (%.4f%%)\n', GCI_cl_fine, GCI_cl_fine*100);
fprintf('Asymptotic ratio:      %.4f (should be ~1.0)\n', asymptotic_cl);

fprintf('\n=========== GCI Results: Cl/Cd ===========\n');
fprintf('     Grid1 Cl/Cd:      %.6f\n', Grid1_clcd);
fprintf('     Grid2 Cl/Cd:      %.6f\n', Grid2_clcd);
fprintf('     Grid3 Cl/Cd:      %.6f\n', Grid3_clcd);
fprintf('     Grid4 Cl/Cd:      %.6f\n', Grid4_clcd);
fprintf('     Grid5 Cl/Cd:      %.6f\n', Grid5_clcd);
fprintf('   f_exact Cl/Cd:      %.6f\n', f_exact_clcd);
fprintf(' GCI_supercoarse:      %.6f (%.4f%%)\n', GCI_clcd_supercoarse, GCI_clcd_supercoarse*100);
fprintf('      GCI_coarse:      %.6f (%.4f%%)\n', GCI_clcd_coarse, GCI_clcd_coarse*100);
fprintf('      GCI_medium:      %.6f (%.4f%%)\n', GCI_clcd_medium, GCI_clcd_medium*100);
fprintf('        GCI_fine:      %.6f (%.4f%%)\n', GCI_clcd_fine, GCI_clcd_fine*100);
fprintf('Asymptotic ratio:      %.4f (should be ~1.0)\n', asymptotic_clcd);

%% Convergence plots for spatial grid
figure(1)
subplot(2,1,1)
plot(h_array, [Grid5_cl, Grid4_cl, Grid3_cl, Grid2_cl, Grid1_cl], 'bo-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cl, 'b*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, h_coarse], [f_exact_cl, Grid3_cl], 'b--', 'LineWidth', 1)
hold off
xlabel('Mesh Spacing h = \surd(1/N)')
ylabel('C_L')
title('Lift Coefficient vs Mesh Spacing')
legend('CFD Solutions', 'Richardson Extrapolation (h→0)', 'Extrapolation trend', ...
    'Location', 'best')
%xlim([-h_supercoarse*0.1, h_supercoarse*1.1])
grid on
gca.XAxis.Exponent = 0;

subplot(2,1,2)
plot(h_array, [Grid5_cd, Grid4_cd, Grid3_cd, Grid2_cd, Grid1_cd], 'ro-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cd, 'r*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, h_coarse], [f_exact_cd, Grid3_cd], 'r--', 'LineWidth', 1)
hold off
xlabel('Mesh Spacing h = \surd(1/N)')
ylabel('C_D')
title('Drag Coefficient vs Mesh Spacing')
legend('CFD Solutions', 'Richardson Extrapolation (h→0)', 'Extrapolation trend', ...
    'Location', 'best')
%xlim([-h_supercoarse*0.1, h_supercoarse*1.1])
grid on
gca.XAxis.Exponent = 0;

%% Convergence plots for temporal grid (raw dt)
dt_array = [dt_fine, dt_medium, dt_coarse, dt_mediumcoarse, dt_supercoarse];

figure(2)
subplot(2,1,1)
plot(dt_array, [Grid5_cl, Grid4_cl, Grid3_cl, Grid2_cl, Grid1_cl], 'bo-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
ax1 = gca;
ax1.XAxis.Exponent = 0;
hold on
plot(0, f_exact_cl, 'b*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, dt_coarse], [f_exact_cl, Grid3_cl], 'b--', 'LineWidth', 1)
hold off
xlabel('\Deltat (seconds)')
ylabel('C_L')
title('Lift Coefficient vs Time Step')
legend('CFD Solutions', 'Richardson Extrapolation (\Deltat→0)', 'Extrapolation trend', ...
    'Location', 'best')
%xlim([-dt_supercoarse*0.1, dt_supercoarse*1.1])
grid on
gca.XAxis.Exponent = 0;

subplot(2,1,2)
plot(dt_array, [Grid5_cd, Grid4_cd, Grid3_cd, Grid2_cd, Grid1_cd], 'ro-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cd, 'r*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, dt_coarse], [f_exact_cd, Grid3_cd], 'r--', 'LineWidth', 1)
hold off
xlabel('\Deltat (seconds)')
ylabel('C_D')
title('Drag Coefficient vs Time Step')
legend('CFD Solutions', 'Richardson Extrapolation (\Deltat→0)', 'Extrapolation trend', ...
    'Location', 'best')
%xlim([-dt_supercoarse*0.1, dt_supercoarse*1.1])
grid on
gca.XAxis.Exponent = 0;

%% Convergence plots for normalized spatial grid
figure(3)
subplot(2,1,1)
plot(h_norm, [Grid5_cl, Grid4_cl, Grid3_cl, Grid2_cl, Grid1_cl], 'bo-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cl, 'b*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, h_norm_coarse], [f_exact_cl, Grid3_cl], 'b--', 'LineWidth', 1)
hold off
xlabel('Normalised Mesh Spacing h/h_{fine}')
ylabel('C_L')
title('Lift Coefficient vs Normalised Mesh Spacing')
legend('CFD Solutions', 'Richardson Extrapolation (h→0)', 'Extrapolation trend', ...
    'Location', 'best')
xlim([-h_norm_supercoarse*0.1, h_norm_supercoarse*1.1])
grid on


subplot(2,1,2)
plot(h_norm, [Grid5_cd, Grid4_cd, Grid3_cd, Grid2_cd, Grid1_cd], 'ro-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cd, 'r*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, h_norm_coarse], [f_exact_cd, Grid3_cd], 'r--', 'LineWidth', 1)
hold off
xlabel('Normalised Mesh Spacing h/h_{fine}')
ylabel('C_D')
title('Drag Coefficient vs Normalised Mesh Spacing')
legend('CFD Solutions', 'Richardson Extrapolation (h→0)', 'Extrapolation trend', ...
    'Location', 'best')
xlim([-h_norm_supercoarse*0.1, h_norm_supercoarse*1.1])
grid on


%% Convergence plots for normalized temporal grid
dt_norm_fine         = dt_fine         / dt_fine;
dt_norm_medium       = dt_medium       / dt_fine;
dt_norm_coarse       = dt_coarse       / dt_fine;
dt_norm_mediumcoarse = dt_mediumcoarse / dt_fine;
dt_norm_supercoarse  = dt_supercoarse  / dt_fine;
dt_norm = [dt_norm_fine, dt_norm_medium, dt_norm_coarse, ...
           dt_norm_mediumcoarse, dt_norm_supercoarse];

figure(4)
subplot(2,1,1)
plot(dt_norm, [Grid5_cl, Grid4_cl, Grid3_cl, Grid2_cl, Grid1_cl], 'bo-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cl, 'b*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, dt_norm_coarse], [f_exact_cl, Grid3_cl], 'b--', 'LineWidth', 1)
hold off
xlabel('Normalised Time Step \Deltat/\Deltat_{fine}')
ylabel('C_L')
title('Lift Coefficient vs Normalised Time Step')
legend('CFD Solutions', 'Richardson Extrapolation (\Deltat→0)', 'Extrapolation trend', ...
    'Location', 'best')
xlim([-dt_norm_supercoarse*0.1, dt_norm_supercoarse*1.1])
grid on

subplot(2,1,2)
plot(dt_norm, [Grid5_cd, Grid4_cd, Grid3_cd, Grid2_cd, Grid1_cd], 'ro-', ...
    'LineWidth', 1.5, 'MarkerSize', 8)
hold on
plot(0, f_exact_cd, 'r*', 'MarkerSize', 12, 'LineWidth', 2)
plot([0, dt_norm_coarse], [f_exact_cd, Grid3_cd], 'r--', 'LineWidth', 1)
hold off
xlabel('Normalised Time Step \Deltat/\Deltat_{fine}')
ylabel('C_D')
title('Drag Coefficient vs Normalised Time Step')
legend('CFD Solutions', 'Richardson Extrapolation (\Deltat→0)', 'Extrapolation trend', ...
    'Location', 'best')
xlim([-dt_norm_supercoarse*0.1, dt_norm_supercoarse*1.1])
grid on